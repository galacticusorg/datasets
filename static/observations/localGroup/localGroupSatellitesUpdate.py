#!/usr/bin/env python3
"""Update bibliographic references in localGroupSatellites.xml via NASA ADS API."""

import os
import re
import sys
import shutil
import unicodedata
import urllib.parse

import requests


def resolve_url(url):
    """Follow redirects and return the final URL, stripping any /abstract suffix."""
    resp = requests.get(url, allow_redirects=True)
    final_url = re.sub(r"/abstract$", "", resp.url)
    return final_url


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: localGroupSatellitesUpdate.py <apiToken>")

    api_token = sys.argv[1]
    known_updated_urls = {}
    bib_codes = {}

    # Stage 1: Update URLs and collect bibcodes.
    with open("localGroupSatellites.xml", "r") as fin, \
         open("localGroupSatellites.xml.stage1", "w") as fout:
        for line in fin:
            # Update URLs pointing to an arXiv paper on NASA ADS.
            m = re.search(r'"(https://ui\.adsabs\.harvard\.edu/abs/\d+arXiv.*?)"', line)
            if m:
                old_url = m.group(1)
                if old_url not in known_updated_urls:
                    new_url = resolve_url(old_url)
                    known_updated_urls[old_url] = new_url
                    print(f"Updating URL '{old_url}' to '{new_url}'")
                line = line.replace(old_url, known_updated_urls[old_url])

            # Update URLs pointing directly to arXiv.
            m = re.search(r'"(https?://arxiv\.org/abs/(\d+\.\d+))"', line)
            if m:
                original_url = m.group(1)
                arxiv_id = m.group(2)
                ads_url = f"https://ui.adsabs.harvard.edu/abs/arXiv:{arxiv_id}"
                if ads_url not in known_updated_urls:
                    new_url = resolve_url(ads_url)
                    known_updated_urls[ads_url] = new_url
                    print(f"Updating URL '{original_url}' to '{new_url}'")
                line = line.replace(original_url, known_updated_urls[ads_url])

            # Collect bibcodes from ADS URLs.
            m = re.search(r'"https://ui\.adsabs\.harvard\.edu/abs/(.*?)"', line)
            if m:
                bib_code = urllib.parse.unquote(m.group(1))
                bib_code = re.sub(r"/abstract$", "", bib_code)
                bib_codes[bib_code] = "unknown"

            fout.write(line)

    # Query NASA ADS bigquery endpoint for bibliographic metadata.
    count = len(bib_codes)
    url = (
        "https://api.adsabs.harvard.edu/v1/search/bigquery"
        f"?q=*:*&rows={count}&fl=bibcode,alternate_bibcode,title,author,year,pub,volume,page"
    )
    headers = {
        "Authorization": f"Bearer {api_token}",
        "Content-Type": "big-query/csv",
    }
    post_body = "bibcode\n" + "\n".join(sorted(bib_codes))
    resp = requests.post(url, headers=headers, data=post_body)
    if resp.status_code != 200:
        sys.exit(f"Failed to retrieve records: {resp.status_code}\n{resp.text}")
    records = resp.json()

    # Resolve alternate bibcodes: if we used an alternate, switch to it.
    for record in records["response"]["docs"]:
        record["canonical_bibcode"] = record["bibcode"]
        for alt in record.get("alternate_bibcode", []):
            if alt in bib_codes:
                record["bibcode"] = alt

    # Journal abbreviation table.
    journal_abbr = {
        "The Astronomical Journal":                                "AJ",
        "The Astrophysical Journal":                               "ApJ",
        "Astronomy and Astrophysics":                              "AA",
        "Acta Astronomica":                                        "Acta Astron.",
        "Monthly Notices of the Royal Astronomical Society":       "MNRAS",
        "arXiv e-prints":                                          "arXiv",
        "Research Notes of the American Astronomical Society":     "RNAAS",
        "Publications of the Astronomical Society of the Pacific": "PASJ",
        "Publications of the Astronomical Society of Japan":       "PASP",
        "Nature":                                                  "Nature",
    }

    bib_codes_canonical = {}
    for record in records["response"]["docs"]:
        pub = record.get("pub", "")
        if pub not in journal_abbr:
            sys.exit(f"No journal abbreviation found for '{pub}'")

        raw_authors = [re.sub(r"([^,]+),.*", r"\1", a) for a in record.get("author", [])]
        n = len(raw_authors)
        if n == 1:
            author = raw_authors[0]
        elif n == 2:
            author = f"{raw_authors[0]} &amp; {raw_authors[1]}"
        elif n == 3:
            author = f"{raw_authors[0]}, {raw_authors[1]} &amp; {raw_authors[2]}"
        else:
            author = f"{raw_authors[0]}, {raw_authors[1]}, {raw_authors[2]} et al."
        author += "; "

        # Strip diacritics via NFKD normalization.
        author = unicodedata.normalize("NFKD", author)
        author = "".join(c for c in author if unicodedata.category(c) != "Mn")

        year    = record["year"] + "; "
        journal = journal_abbr[pub] + "; "
        volume  = record["volume"] + "; " if "volume" in record else ""
        page    = record["page"][0]

        bib_codes[record["bibcode"]] = author + year + journal + volume + page
        bib_codes_canonical[record["bibcode"]] = record["canonical_bibcode"]

    # Stage 2: Update reference attributes in the XML.
    with open("localGroupSatellites.xml.stage1", "r") as fin, \
         open("localGroupSatellites.xml.stage2", "w") as fout:
        for line in fin:
            m = re.search(r'"https://ui\.adsabs\.harvard\.edu/abs/(.*)"', line)
            if m:
                bib_code = m.group(1)
                if bib_code in bib_codes and bib_code in bib_codes_canonical:
                    canonical = bib_codes_canonical[bib_code]
                    line = re.sub(r'\sreference="[^"]+"',
                                  f' reference="{bib_codes[bib_code]}"', line)
                    line = re.sub(r'\sreferenceURL="[^"]+"',
                                  f' referenceURL="https://ui.adsabs.harvard.edu/abs/{canonical}"',
                                  line)
            fout.write(line)

    shutil.move("localGroupSatellites.xml.stage2", "localGroupSatellites.xml")
    os.unlink("localGroupSatellites.xml.stage1")


if __name__ == "__main__":
    main()
