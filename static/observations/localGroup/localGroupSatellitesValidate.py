#!/usr/bin/env python3
"""Validate localGroupSatellites.xml against its XSD schema, and check that the positions it records are consistent.

Most galaxies record their position in several forms, some of which are derived from others by localGroupDBUpdate in
Galacticus (source/interface/Local_Group_database.F90). That function does not recompute a derived field which is already
present, so a correction to one form does not reach the others unless those fields are removed and regenerated. The checks
here are limited to what can be established from this file alone: they find disagreements, not values which are wrong in
every form.
"""

import math
import sys
import xml.etree.ElementTree as ET

import xmlschema

DATABASE = "localGroupSatellites.xml"
SCHEMA = "localGroupSatellites.xsd"

# J2000 position of the north Galactic pole, and the Galactic longitude of the north celestial pole, in degrees.
POLE_RIGHT_ASCENSION = 192.85948
POLE_DECLINATION = 27.12825
POLE_LONGITUDE = 122.93192

# Largest allowed separation between (l,b) and (RA,Dec), in degrees. Stored Galactic coordinates are often given to 0.1
# degrees: the largest separation among current entries is 0.07 degrees.
TOLERANCE_GALACTIC = 0.1
# Largest allowed difference between the sexagesimal and decimal forms of a coordinate, in degrees.
TOLERANCE_SEXAGESIMAL = 1.0e-4
# Largest allowed fractional difference for fields derived by localGroupDBUpdate, which are written to 8 significant figures,
# and the corresponding absolute allowance in Mpc, for distances close to zero.
TOLERANCE_DERIVED_FRACTIONAL = 1.0e-6
TOLERANCE_DERIVED_ABSOLUTE = 1.0e-7
# Largest allowed fractional difference between the heliocentric distance and that implied by the distance modulus. These can
# be separate measurements, so this is intended only to find gross errors, such as those in the units or the decimal point.
TOLERANCE_DISTANCE_MODULUS = 0.15

REGENERATE = "remove the derived fields and regenerate them with localGroupDBUpdate"


def value(galaxy, tag):
    """Return the value of the named element as a float, or None if the element is absent."""
    element = galaxy.find(tag)
    return None if element is None else float(element.get("value"))


def vector(galaxy, tag):
    """Return the value of the named element as a list of floats, or None if the element is absent."""
    element = galaxy.find(tag)
    return None if element is None else [float(x) for x in element.get("value").split()]


def right_ascension_sexagesimal(galaxy):
    """Return the sexagesimal right ascension in degrees, or None if it is absent."""
    element = galaxy.find("rightAscensionSexagesimal")
    if element is None:
        return None
    return 15.0 * (float(element.get("hours")) + float(element.get("minutes")) / 60.0 + float(element.get("seconds")) / 3600.0)


def declination_sexagesimal(galaxy):
    """Return the sexagesimal declination in degrees, or None if it is absent. The sign is taken from the text, as a
    declination between 0 and -1 degrees is written with degrees of "-0"."""
    element = galaxy.find("declinationSexagesimal")
    if element is None:
        return None
    degrees = element.get("degrees").strip()
    magnitude = abs(float(degrees)) + float(element.get("arcminutes")) / 60.0 + float(element.get("arcseconds")) / 3600.0
    return -magnitude if degrees.startswith("-") else magnitude


def equatorial_to_galactic(right_ascension, declination):
    """Convert J2000 equatorial coordinates to Galactic longitude and latitude, all in degrees."""
    alpha = math.radians(right_ascension - POLE_RIGHT_ASCENSION)
    delta = math.radians(declination)
    pole = math.radians(POLE_DECLINATION)
    latitude = math.asin(math.sin(pole) * math.sin(delta) + math.cos(pole) * math.cos(delta) * math.cos(alpha))
    longitude = POLE_LONGITUDE - math.degrees(
        math.atan2(math.cos(delta) * math.sin(alpha), math.cos(pole) * math.sin(delta) - math.sin(pole) * math.cos(delta) * math.cos(alpha))
    )
    return longitude % 360.0, math.degrees(latitude)


def separation(longitude1, latitude1, longitude2, latitude2):
    """Return the angle between two directions on the sky, all in degrees."""
    phi1, theta1, phi2, theta2 = map(math.radians, (longitude1, latitude1, longitude2, latitude2))
    delta = phi2 - phi1
    numerator = math.hypot(math.cos(theta2) * math.sin(delta), math.cos(theta1) * math.sin(theta2) - math.sin(theta1) * math.cos(theta2) * math.cos(delta))
    return math.degrees(math.atan2(numerator, math.sin(theta1) * math.sin(theta2) + math.cos(theta1) * math.cos(theta2) * math.cos(delta)))


def check_ranges(galaxy):
    """Check that coordinates lie within their allowed ranges."""
    problems = []
    for tag, low, high in (
        ("rightAscensionDecimal", 0.0, 360.0),
        ("declinationDecimal", -90.0, 90.0),
        ("longitudeGalactic", 0.0, 360.0),
        ("latitudeGalactic", -90.0, 90.0),
    ):
        coordinate = value(galaxy, tag)
        if coordinate is not None and not low <= coordinate <= high:
            problems.append(f"{tag} of {coordinate} is outside the range {low} to {high}")
    return problems


def check_sexagesimal(galaxy):
    """Check that the sexagesimal and decimal forms of each coordinate agree."""
    problems = []
    for tag, sexagesimal in (
        ("rightAscensionDecimal", right_ascension_sexagesimal(galaxy)),
        ("declinationDecimal", declination_sexagesimal(galaxy)),
    ):
        decimal = value(galaxy, tag)
        if decimal is not None and sexagesimal is not None and abs(decimal - sexagesimal) > TOLERANCE_SEXAGESIMAL:
            problems.append(f"{tag} of {decimal} disagrees with the sexagesimal value of {sexagesimal:.6f}")
    return problems


def check_galactic(galaxy):
    """Check that the Galactic coordinates agree with the right ascension and declination."""
    right_ascension, declination = value(galaxy, "rightAscensionDecimal"), value(galaxy, "declinationDecimal")
    longitude, latitude = value(galaxy, "longitudeGalactic"), value(galaxy, "latitudeGalactic")
    if None in (right_ascension, declination, longitude, latitude):
        return []
    expected = equatorial_to_galactic(right_ascension, declination)
    angle = separation(longitude, latitude, *expected)
    if angle <= TOLERANCE_GALACTIC:
        return []
    return [f"(l,b) = ({longitude},{latitude}) is {angle:.2f} degrees from ({expected[0]:.2f},{expected[1]:.2f}), computed from (RA,Dec) = ({right_ascension},{declination})"]


def differs(first, second):
    """Return true if two distances, in Mpc, differ by more than allowed for derived fields."""
    return abs(first - second) > TOLERANCE_DERIVED_FRACTIONAL * max(abs(first), abs(second)) + TOLERANCE_DERIVED_ABSOLUTE


def check_cartesian(galaxy):
    """Check that the Cartesian position is that computed from the right ascension, declination, and distance."""
    element = galaxy.find("positionHeliocentricCartesian")
    right_ascension, declination = value(galaxy, "rightAscensionDecimal"), value(galaxy, "declinationDecimal")
    distance = value(galaxy, "distanceHeliocentric")
    if element is None or None in (right_ascension, declination, distance):
        return []
    alpha, delta = math.radians(right_ascension), math.radians(declination)
    expected = [distance * math.cos(delta) * math.cos(alpha), distance * math.cos(delta) * math.sin(alpha), distance * math.sin(delta)]
    stored = vector(galaxy, "positionHeliocentricCartesian")
    if any(differs(x, y) for x, y in zip(stored, expected)):
        return [f"positionHeliocentricCartesian of {stored} disagrees with {[float(f'{x:.8e}') for x in expected]}, computed from (RA,Dec) and distanceHeliocentric: {REGENERATE}"]
    return []


def check_distances(galaxy, positions):
    """Check that the distances to the Milky Way and M31 are those computed from the Cartesian positions. A distance which
    has its own reference may be a published value, and so is not checked."""
    problems = []
    position = vector(galaxy, "positionHeliocentricCartesian")
    if position is None:
        return problems
    for tag, reference in positions.items():
        element = galaxy.find(tag)
        if element is None or element.get("reference") is not None:
            continue
        stored, expected = float(element.get("value")), math.dist(position, reference)
        if differs(stored, expected):
            problems.append(f"{tag} of {stored} disagrees with {expected:.8e}, computed from the Cartesian positions: {REGENERATE}")
    return problems


def check_distance_modulus(galaxy):
    """Check that the heliocentric distance is not grossly inconsistent with the distance modulus."""
    distance, modulus = value(galaxy, "distanceHeliocentric"), value(galaxy, "distanceModulus")
    if None in (distance, modulus):
        return []
    expected = 10.0 ** ((modulus - 25.0) / 5.0)
    if abs(distance / expected - 1.0) > TOLERANCE_DISTANCE_MODULUS:
        return [f"distanceHeliocentric of {distance} Mpc disagrees with {expected:.6g} Mpc, computed from distanceModulus of {modulus}"]
    return []


def find_problems(root):
    """Return a list of descriptions of inconsistent positions."""
    galaxies = list(root.iter("galaxy"))
    names = {galaxy.get("name"): galaxy for galaxy in galaxies}
    positions = {
        "distanceMilkyWay": vector(names["The Galaxy"], "positionHeliocentricCartesian"),
        "distanceM31": vector(names["Andromeda"], "positionHeliocentricCartesian"),
    }
    problems = []
    for galaxy in galaxies:
        for problem in (
            check_ranges(galaxy) + check_sexagesimal(galaxy) + check_galactic(galaxy) + check_cartesian(galaxy) + check_distances(galaxy, positions) + check_distance_modulus(galaxy)
        ):
            problems.append(f"{galaxy.get('name')}: {problem}")
    return problems


def main():
    schema = xmlschema.XMLSchema11(SCHEMA)
    schema.validate(DATABASE)
    root = ET.parse(DATABASE).getroot()
    problems = find_problems(root)
    for problem in problems:
        print(problem)
    print(f"Checked positions of {sum(1 for _ in root.iter('galaxy'))} galaxies: {len(problems)} problems found.")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
