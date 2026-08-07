#!/usr/bin/env python3
"""Record membership of the DELVE Milky Way Census I in the Local Group satellite database.

The DELVE Milky Way Census I (Tan, Drlica-Wagner et al. 2025; arXiv:2509.12313) defines a sample of
Milky Way satellites recovered above a uniform detection threshold across the combined footprints of
DES Y6, DELVE DR3, and Pan-STARRS1 DR1, for which a well-defined observational selection function is
available. That sample -- not the set of galaxies discovered by any particular survey -- is the
appropriate target for comparison with models to which that selection function has been applied.

This script reads the census list released with that paper and adds a `<census>` element to each
matching galaxy in `localGroupSatellites.xml`, recording membership and the census classification
("dwarf" for confirmed dwarf galaxies, "probableDwarf" for probable ones). It is idempotent: any
existing `<census>` elements for this census are replaced.

Where a census member has no `<magnitudeAbsoluteV>` in the database, the value compiled by the
census is also added, attributed to the source from which the census took it. Existing values are
never overwritten -- the database's own compilation takes precedence.

The census list is `data/offical_delve_census_refgalaxy_list.csv` of

   https://github.com/delve-survey/delve_mw_census

Usage:

   ./localGroupSatellitesCensusUpdate.py offical_delve_census_refgalaxy_list.csv

Andrew Benson <abenson@carnegiescience.edu>
"""

import csv
import re
import sys
import xml.etree.ElementTree as ET

# Name of this census, as recorded in the database.
censusName = "DELVE-MW-I"

reference = "Tan, Drlica-Wagner et al.; 2025; arXiv:2509.12313"
referenceURL = "https://arxiv.org/abs/2509.12313"

# Census classifications, and the corresponding names used in the database.
classifications = {"D": "dwarf", "PD": "probableDwarf"}

# Census galaxies whose names in the census list can not be matched to a database name or alias by
# normalization alone.
synonyms = {
    "Pictor I": "Pictoris I",
    "Pictor II": "Pictoris II",
    "Sagittarius": "Sagittarius dSph",
    "Segue 1": "Segue (I)",
    "Segue 2": "Segue II",
    "Sextans": "Sextans (I)",
}

fileNameDatabase = "localGroupSatellites.xml"


def normalize(name):
    """Reduce a galaxy name to a form suitable for matching."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def citation(reference):
    """Split a census reference into a human-readable citation and an ADS URL.

    References in the census list are an author name concatenated with an ADS bibcode, for example
    `Torrealba2018MNRAS.475.5085T`.
    """
    match = re.match(r"^(?P<author>[^\d]+)(?P<bibcode>\d{4}.*)$", reference)
    if match is None:
        sys.exit(f"can not parse census reference '{reference}'")
    author = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", match.group("author"))
    return (
        f"{author} ({match.group('bibcode')[0:4]})",
        f"https://ui.adsabs.harvard.edu/abs/{match.group('bibcode')}",
    )


def main(fileNameCensus):
    # Read the census, retaining only those systems which are members.
    census = {
        row["name"]: row
        for row in csv.DictReader(open(fileNameCensus))
        if row["delve_census"] == "True"
    }
    unknown = set(row["dwarf_class"] for row in census.values()) - set(classifications.keys())
    if unknown:
        sys.exit(f"unknown census classification(s): {', '.join(sorted(unknown))}")

    # Build a map from all names and aliases in the database to the galaxy name.
    names = {}
    for galaxy in ET.parse(fileNameDatabase).getroot().find("galaxies").findall("galaxy"):
        for name in [galaxy.get("name")] + [alias.get("value") for alias in galaxy.findall("alias")]:
            if name is not None:
                names.setdefault(normalize(name), galaxy.get("name"))

    # Match each census galaxy to a galaxy in the database.
    members = {}
    magnitudesCensus = {}
    for name, row in census.items():
        nameNormalized = normalize(synonyms.get(name, name))
        if nameNormalized not in names:
            sys.exit(f"census galaxy '{name}' has no match in '{fileNameDatabase}'")
        members[names[nameNormalized]] = row["dwarf_class"]
        magnitudesCensus[names[nameNormalized]] = (
            float(row["M_V"]),
            row["M_V_err"],
        ) + citation(row["ref_m_v"])
    print(f"Matched {len(members)} of {len(census)} census galaxies to the database.")

    # Insert the census element into each member. The file is edited as text, rather than being
    # re-serialized from a parse tree, so that the formatting of unmodified entries is preserved.
    # Determine which members have no absolute magnitude in the database, and so should take the
    # value compiled by the census.
    magnitudes = {}
    for galaxy in ET.parse(fileNameDatabase).getroot().find("galaxies").findall("galaxy"):
        if galaxy.get("name") in members and galaxy.find("magnitudeAbsoluteV") is None:
            magnitudes[galaxy.get("name")] = magnitudesCensus[galaxy.get("name")]
    if magnitudes:
        print(f"Adding absolute magnitudes compiled by the census for: {', '.join(sorted(magnitudes))}.")

    database = open(fileNameDatabase).read()
    # Remove any existing entries for this census, so that the script may be re-run.
    database = re.sub(rf'^[ \t]*<census value="{censusName}"[^>]*/>\n', "", database, flags=re.MULTILINE)
    countInserted = 0

    def insert(match):
        nonlocal countInserted
        name = match.group("name")
        if name not in members:
            return match.group(0)
        countInserted += 1
        indent = match.group("indent")
        element = (
            f'{indent}<census value="{censusName}"'
            f' classification="{classifications[members[name]]}"'
            f' reference="{reference}" referenceURL="{referenceURL}" />\n'
        )
        if name in magnitudes:
            magnitude, uncertainty, referenceMagnitude, referenceMagnitudeURL = magnitudes[name]
            element += (
                f'{indent}<magnitudeAbsoluteV value="{magnitude:+.3f}" uncertainty="{uncertainty}"'
                f' reference="{referenceMagnitude}" referenceURL="{referenceMagnitudeURL}"/>\n'
            )
        return match.group(0) + element

    # Match the opening tag of each galaxy, together with the indentation of the line which follows
    # it, so that the new element can be inserted as the first child at the correct indentation.
    database, count = re.subn(
        r'(?P<open><galaxy name="(?P<name>[^"]*)"[^>]*>\n)(?=(?P<indent>[ \t]*)<)',
        insert,
        database,
    )
    if countInserted != len(members):
        sys.exit(f"inserted {countInserted} census elements, but matched {len(members)} galaxies")
    open(fileNameDatabase, "w").write(database)
    print(f"Wrote {countInserted} census entries to '{fileNameDatabase}'.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
