#!/usr/bin/env python3
"""Update bibliographic references in localGroupSatellites.xml via NASA ADS API."""

import os
import re
import sys
import shutil
import unicodedata
import urllib.parse

import requests


# Matches an ADS abstract URL, capturing the bibcode. Both the current host and the
# older 'adsabs.harvard.edu' form are accepted, as both appear in the database.
ads_url_pattern = r'"(https?://(?:ui\.)?adsabs\.harvard\.edu/abs/([^"]+))"'


def resolve_url(url):
    """Follow redirects and return the final URL, stripping any /abstract suffix."""
    resp = requests.get(url, allow_redirects=True)
    final_url = re.sub(r"/abstract$", "", resp.url)
    return final_url


def decode_bib_code(text):
    """Extract a bibcode from the path of an ADS URL.

    Bibcodes containing an ampersand (e.g. '2019A&A...625A...2K') are percent-encoded
    in the URL, and some entries in the database are encoded more than once, so decode
    repeatedly until the result stops changing.
    """
    bib_code = re.sub(r"/abstract$", "", text)
    while True:
        decoded = urllib.parse.unquote(bib_code)
        if decoded == bib_code:
            return bib_code
        bib_code = decoded


# Journal abbreviation table, keyed on the exact ADS "pub" field.
journal_abbr = {
    "The Astronomical Journal":                                "AJ",
    "The Astrophysical Journal":                               "ApJ",
    "The Astrophysical Journal Letters":                       "ApJL",
    "The Astrophysical Journal Supplement Series":             "ApJS",
    "Astronomy and Astrophysics":                              "AA",
    "Astronomy and Astrophysics Supplement Series":            "AAS",
    "Annual Review of Astronomy and Astrophysics":             "ARAA",
    "Acta Astronomica":                                        "Acta Astron.",
    "Monthly Notices of the Royal Astronomical Society":       "MNRAS",
    "arXiv e-prints":                                          "arXiv",
    "Research Notes of the American Astronomical Society":     "RNAAS",
    "Publications of the Astronomical Society of the Pacific": "PASP",
    "Publications of the Astronomical Society of Japan":       "PASJ",
    "Reviews of Modern Physics":                               "RvMP",
    "Nature":                                                  "Nature",
    "Nature Astronomy":                                        "Nature Astron.",
}

# Books and catalogs, for which the ADS "pub" field holds the full title. Matched by
# regex as these titles are long and can include volume-by-volume detail.
publication_abbr = [
    (r"^Third Reference Catalogue of Bright Galaxies", "RC3"),
]


def abbreviate(pub):
    """Return the abbreviation for a publication name, or None if it is unknown."""
    if pub in journal_abbr:
        return journal_abbr[pub]
    for pattern, abbr in publication_abbr:
        if re.search(pattern, pub):
            return abbr
    return None


def format_reference(record, journal):
    """Build the reference string for an ADS record: authors; year; journal; volume; page."""
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

    # Strip diacritics via NFKD normalization.
    author = unicodedata.normalize("NFKD", author)
    author = "".join(c for c in author if unicodedata.category(c) != "Mn")

    # Books and catalogs carry no volume or page, so include only those parts present.
    parts = [author, record["year"], journal]
    if record.get("volume"):
        parts.append(record["volume"])
    if record.get("page"):
        parts.append(record["page"][0])
    return "; ".join(parts)


def update_line(line, bib_codes, bib_codes_canonical):
    """Rewrite the reference and referenceURL attributes on a line citing an ADS record."""
    m = re.search(ads_url_pattern, line)
    if not m:
        return line

    bib_code = decode_bib_code(m.group(2))
    if bib_code not in bib_codes or bib_code not in bib_codes_canonical:
        return line

    # Percent-encode the bibcode, as an ampersand can not appear literally in the URL.
    # The colon of an 'arXiv:' bibcode is left as is, matching the URLs that ADS serves.
    canonical = urllib.parse.quote(bib_codes_canonical[bib_code], safe=":")
    line = re.sub(r'\sreference="[^"]+"',
                  f' reference="{bib_codes[bib_code]}"', line)
    line = re.sub(r'\sreferenceURL="[^"]+"',
                  f' referenceURL="https://ui.adsabs.harvard.edu/abs/{canonical}"', line)
    return line


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
            m = re.search(r'"(https?://(?:ui\.)?adsabs\.harvard\.edu/abs/\d+arXiv[^"]*?)"', line)
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
            m = re.search(ads_url_pattern, line)
            if m:
                bib_codes[decode_bib_code(m.group(2))] = "unknown"

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

    # Accumulate any unknown publications so that a single run reports all of them,
    # rather than failing on the first and hiding the rest until the next run.
    bib_codes_canonical = {}
    unknown_pubs = {}
    for record in records["response"]["docs"]:
        pub = record.get("pub", "")
        journal = abbreviate(pub)
        if journal is None:
            unknown_pubs.setdefault(pub, []).append(record["bibcode"])
            continue

        bib_codes[record["bibcode"]] = format_reference(record, journal)
        bib_codes_canonical[record["bibcode"]] = record["canonical_bibcode"]

    if unknown_pubs:
        report = "\n".join(f"  '{pub}' ({', '.join(sorted(codes))})"
                           for pub, codes in sorted(unknown_pubs.items()))
        sys.exit(f"No journal abbreviation found for:\n{report}")

    # Stage 2: Update reference attributes in the XML.
    with open("localGroupSatellites.xml.stage1", "r") as fin, \
         open("localGroupSatellites.xml.stage2", "w") as fout:
        for line in fin:
            fout.write(update_line(line, bib_codes, bib_codes_canonical))

    shutil.move("localGroupSatellites.xml.stage2", "localGroupSatellites.xml")
    os.unlink("localGroupSatellites.xml.stage1")


if __name__ == "__main__":
    main()
