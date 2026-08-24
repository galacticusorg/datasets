#!/usr/bin/env python3
"""Convert the Ferrara et al. (1999) dust attenuation atlas from XML to HDF5.

The atlas is distributed here as XML, one file per combination of dust type and dust-to-stellar scale height ratio.
The dust compendium of Benson (2018) tabulates the same kind of quantity -- the fraction of light escaping a galaxy
as a function of wavelength, inclination, optical depth, and spheroid size -- and does so in HDF5, with the same four
axes in the same order. Converting to that layout lets one reader serve both atlases, and moves a four-deep nested
XML parse out of the model code.

The layout written here follows the compendium:

    /wavelength           (nWavelength)                                             Angstroms
    /inclination          (nInclination)                                            degrees
    /opticalDepth         (nOpticalDepth)                                           V band, face-on, through center
    /spheroidScaleRadial  (nSpheroidScaleRadial)                                    units of the disk scale length
    /attenuationDisk      (nWavelength, nInclination, nOpticalDepth)                transmitted fraction
    /attenuationSpheroid  (nWavelength, nInclination, nOpticalDepth, nSpheroidScaleRadial)

Note on `/spheroidScaleRadial`. The compendium's spheroids are Hernquist profiles and this axis is their scale
radius; Ferrara et al. use Jaffe profiles and tabulate against the spheroid effective radius, eta = r_e/r_star, with
r_star the disk radial scale length. For a Jaffe profile rho(r) = (rho_0/4pi) (r/r_0)^-2 (1+r/r_0)^-2 the enclosed
mass is M(r)/M = (r/r_0)/(1+r/r_0), so the half-mass radius is exactly r_0 and eta is the Jaffe scale radius in units
of the disk scale length. The two atlases therefore share this axis's *format* but not its profile: a model galaxy
must be converted onto it through a radius which means the same thing in both, which is the half-mass radius. That
conversion belongs in the model, not here, and the profile assumed is recorded in the `spheroidProfile` attribute so
that it cannot be lost.

The conversion is faithful: no points are added, removed, or altered. In particular the tabulation begins at an
optical depth of 0.1, and although Ferrara et al. state a grid including tau_V = 0, the distributed files do not
contain it. Deciding what to do below the tabulated range is left to the model, where it is visible.

Andrew Benson, with assistance from Claude.
"""
import os
import sys

import h5py
import numpy as np
from lxml import etree

# Dust types and scale height ratios, as they appear in the file names.
dustTypes           = ("MilkyWay", "SmallMagellanicCloud")
scaleHeightRatios   = ("0.4", "1.0", "2.5")


def readAtlas(fileName):
    """Read one atlas file, returning its axes and the two attenuation tables."""
    root = etree.parse(fileName).getroot()

    wavelength = np.array([float(element.text) for element in root.find("wavelengths")])

    components = {component.findtext("name"): component for component in root.findall("components")}
    for name in ("disk", "bulge"):
        if name not in components:
            raise ValueError(f"{fileName}: no '{name}' component")

    def axes(inclinationElements):
        inclination  = np.array([float(element.findtext("angle")) for element in inclinationElements])
        opticalDepth = np.array([float(element.findtext("tau"  )) for element in inclinationElements[0].findall("opticalDepth")])
        return inclination, opticalDepth

    def table(inclinationElements, inclination, opticalDepth):
        values = np.zeros((wavelength.size, inclination.size, opticalDepth.size))
        for i, inclinationElement in enumerate(inclinationElements):
            depthElements = inclinationElement.findall("opticalDepth")
            if len(depthElements) != opticalDepth.size:
                raise ValueError(f"{fileName}: ragged optical depth axis")
            for j, depthElement in enumerate(depthElements):
                if float(depthElement.findtext("tau")) != opticalDepth[j]:
                    raise ValueError(f"{fileName}: optical depth axis is not shared between inclinations")
                attenuations = depthElement.findall("attenuation")
                if len(attenuations) != wavelength.size:
                    raise ValueError(f"{fileName}: ragged wavelength axis")
                values[:, i, j] = [float(element.text) for element in attenuations]
        return values

    # The disk is tabulated directly against inclination; the bulge against bulge size first.
    inclinationElements                = components["disk"].findall("inclination")
    inclination, opticalDepth          = axes(inclinationElements)
    attenuationDisk                    = table(inclinationElements, inclination, opticalDepth)

    sizeElements                       = components["bulge"].findall("bulgeSize")
    spheroidScaleRadial                = np.array([float(element.findtext("size")) for element in sizeElements])
    attenuationSpheroid                = np.zeros((wavelength.size, inclination.size, opticalDepth.size, spheroidScaleRadial.size))
    for k, sizeElement in enumerate(sizeElements):
        inclinationElementsSpheroid    = sizeElement.findall("inclination")
        inclinationSpheroid, depthSpheroid = axes(inclinationElementsSpheroid)
        if not np.array_equal(inclinationSpheroid, inclination) or not np.array_equal(depthSpheroid, opticalDepth):
            raise ValueError(f"{fileName}: disk and bulge are tabulated on different axes")
        attenuationSpheroid[:, :, :, k] = table(inclinationElementsSpheroid, inclination, opticalDepth)

    return dict(
        wavelength          = wavelength         ,
        inclination         = inclination        ,
        opticalDepth        = opticalDepth       ,
        spheroidScaleRadial = spheroidScaleRadial,
        attenuationDisk     = attenuationDisk    ,
        attenuationSpheroid = attenuationSpheroid,
        description         = (root.findtext("description") or "").strip(),
        scaleHeightStellar  = float(root.findtext("hz_rs" )),
        scaleHeightDust     = float(root.findtext("hzd_hz")),
    )


def main():
    path = os.path.join(os.environ.get("GALACTICUS_DATA_PATH", "."), "static", "dust", "atlasFerrara2000")
    if not os.path.isdir(path):
        sys.exit(f"atlas directory not found: {path}")

    for dustType in dustTypes:
        for scaleHeightRatio in scaleHeightRatios:
            stem     = f"attenuations_{dustType}_dustHeightRatio{scaleHeightRatio}"
            fileXML  = os.path.join(path, stem + ".xml" )
            fileHDF5 = os.path.join(path, stem + ".hdf5")
            if not os.path.exists(fileXML):
                sys.exit(f"missing atlas file: {fileXML}")
            atlas = readAtlas(fileXML)

            # The tabulated values are transmitted fractions. A few percent of them slightly exceed unity -- up to
            # 1.03 -- and that is physical rather than a misread: they occur only at low optical depth and low
            # inclination, where scattering redirects more light into the line of sight than the dust removes from
            # it, so the galaxy is brighter at that angle than it would be with no dust at all. Guard against a
            # genuinely misread file with a bound loose enough to admit that.
            for name in ("attenuationDisk", "attenuationSpheroid"):
                if atlas[name].min() < 0.0 or atlas[name].max() > 1.1:
                    sys.exit(f"{fileXML}: {name} lies outside [0,1.1]: [{atlas[name].min()},{atlas[name].max()}]")

            with h5py.File(fileHDF5, "w") as file:
                for name in ("wavelength", "inclination", "opticalDepth", "spheroidScaleRadial",
                             "attenuationDisk", "attenuationSpheroid"):
                    file.create_dataset(name, data=atlas[name])
                file["wavelength"         ].attrs["units"      ] = "Angstroms"
                file["inclination"        ].attrs["units"      ] = "degrees"
                file["opticalDepth"       ].attrs["description"] = "V-band optical depth, face-on, through the galaxy center"
                file["spheroidScaleRadial"].attrs["description"] = "Jaffe scale radius, equal to the half-mass radius, in units of the disk scale length"
                file.attrs["source"            ] = "Ferrara et al. (1999), ApJS, 123, 437"
                file.attrs["description"       ] = atlas["description"]
                file.attrs["spheroidProfile"   ] = "Jaffe"
                file.attrs["diskProfile"       ] = "exponential"
                file.attrs["scaleHeightStellar"] = atlas["scaleHeightStellar"]
                file.attrs["scaleHeightDust"   ] = atlas["scaleHeightDust"   ]
                file.attrs["convertedFrom"     ] = os.path.basename(fileXML)
            print(f"{stem}: wavelength {atlas['wavelength'].size}, inclination {atlas['inclination'].size}, "
                  f"opticalDepth {atlas['opticalDepth'].size}, spheroidScaleRadial {atlas['spheroidScaleRadial'].size}")


if __name__ == "__main__":
    main()
