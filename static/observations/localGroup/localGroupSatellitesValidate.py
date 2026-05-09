#!/usr/bin/env python3
"""Validate localGroupSatellites.xml against its XSD schema."""

import xmlschema


def main():
    schema = xmlschema.XMLSchema11("localGroupSatellites.xsd")
    schema.validate("localGroupSatellites.xml")


if __name__ == "__main__":
    main()
