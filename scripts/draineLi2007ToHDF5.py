#!/usr/bin/env python3
"""Convert the Draine & Li (2007) dust emission spectra to HDF5.

The spectra are distributed at https://www.astro.princeton.edu/~draine/dust/irem.html as the tarball
irem4/DL07spec.tgz, one ASCII file per grain model and radiation field, named Uumin/Uumin_umax_model.txt. Each gives
the emission per H nucleon of dust heated by starlight with the spectrum of Mathis, Mezger & Panagia (1983) scaled by
U, where the dust mass is distributed over U as

    dM = const U^-2 dU  for umin <= U <= umax,

or, for umax = umin, is all heated by a single intensity. Each file's header gives umin, umax, <U>, and the total
power radiated per H nucleon; its body tabulates nu dP/dnu (erg/s per H) at 1001 wavelengths from 1 cm to 1 micron.

Usage:

    draineLi2007ToHDF5.py DL07spec.tgz

The layout written, to $GALACTICUS_DATA_PATH/static/dust/emission/draineLi2007.hdf5, is:

    /wavelength                 (nWavelength)                                Angstroms, ascending
    /intensityMinimum           (nIntensityMinimum)                          umin
    /intensityMaximum           (nIntensityMaximum)                          umax of the power-law distributions
    /<model>/emissivitySingle   (nIntensityMinimum, nWavelength)             nu dP/dnu, erg/s per H, for U = umin
    /<model>/emissivityPowerLaw (nIntensityMaximum, nIntensityMinimum, nWavelength)
                                                                             nu dP/dnu, erg/s per H, for umin <= U <= umax
    /<model>/powerSingle        (nIntensityMinimum)                          total power, erg/s per H, from the headers
    /<model>/powerPowerLaw      (nIntensityMaximum, nIntensityMinimum)

with the model attributes `massDustPerHydrogen` and `fractionMassPAH` from Draine & Li (2007) Table 3. Groups are named
as the model is enumerated in Galacticus (e.g. `milkyWay60` for MW3.1_60).

The distributed grid differs from the tarball's README, and only the complete, rectangular part of it is written:

  * umin takes 22 values from 0.10 to 25 (the README's 10.0 is absent);
  * umax takes the values 1e3, 1e4, 1e5, 1e6 (the README's 1e7 is absent, and every umax = 1e2 file is an empty stub);
  * single-intensity spectra with U = 1e2 to 3e5 exist only for the MW3.1 models (for LMC2 and SMC they are stubs), so
    they are not written.

Checked on conversion: every file written has the full 1001-point spectrum on a common wavelength grid, with positive
values; its tabulated spectrum integrates to the power in its header to within 0.5%; and the <U> in its header
matches <U> = ln(umax/umin)/(1/umin-1/umax) (or umin) to within 0.5%. The conversion is otherwise faithful: no
points are added, removed, or altered.

Andrew Benson, with assistance from Claude.
"""
import os
import re
import sys
import tarfile

import h5py
import numpy as np

micronsToAngstroms = 1.0e4

# Grain models: the name used in the distributed files, the label used in Galacticus, and, from Draine & Li (2007)
# Table 3, the PAH mass fraction and the dust-to-hydrogen mass ratio.
models = (
    ("MW3.1_00", "milkyWay00"            , 0.0047, 0.0100 ),
    ("MW3.1_10", "milkyWay10"            , 0.0112, 0.0100 ),
    ("MW3.1_20", "milkyWay20"            , 0.0177, 0.0101 ),
    ("MW3.1_30", "milkyWay30"            , 0.0250, 0.0102 ),
    ("MW3.1_40", "milkyWay40"            , 0.0319, 0.0102 ),
    ("MW3.1_50", "milkyWay50"            , 0.0390, 0.0103 ),
    ("MW3.1_60", "milkyWay60"            , 0.0458, 0.0104 ),
    ("LMC2_00" , "largeMagellanicCloud00", 0.0075, 0.00343),
    ("LMC2_05" , "largeMagellanicCloud05", 0.0149, 0.00344),
    ("LMC2_10" , "largeMagellanicCloud10", 0.0237, 0.00359),
    ("smc"     , "smallMagellanicCloud"  , 0.0010, 0.00206),
)
intensityMaximum = np.array([1.0e3, 1.0e4, 1.0e5, 1.0e6])
countWavelength  = 1001
tolerance        = 5.0e-3


def parse(text, name):
    """Parse one spectrum file, returning None for an empty stub."""
    lines = text.splitlines()
    try:
        power         = float(next(line for line in lines if "power/H" in line).split()[0])
        intensityMean = float(next(line for line in lines if "<U>"     in line).split()[0])
        umin, umax    = (float(value) for value in next(line for line in lines if "Umin , Umax" in line).split()[0:2])
        start         = next(i for i, line in enumerate(lines) if line.startswith("lambda"))
    except StopIteration:
        return None
    rows     = [line.split() for line in lines[start+2:]]
    spectrum = np.array([[float(value) for value in row[0:2]] for row in rows if len(row) == 3])
    if spectrum.shape != (countWavelength, 2):
        sys.exit(f"{name}: expected {countWavelength} spectrum rows, found {spectrum.shape[0]}")
    return dict(power=power, intensityMean=intensityMean, umin=umin, umax=umax,
                wavelength=spectrum[::-1, 0]*micronsToAngstroms, emissivity=spectrum[::-1, 1])


def check(record, name):
    """Check a record's spectrum against its own header."""
    if not np.all(record["emissivity"] > 0.0):
        sys.exit(f"{name}: non-positive emissivity")
    frequencyLogarithmic = -np.log(record["wavelength"])
    power                = np.trapezoid(record["emissivity"], frequencyLogarithmic)
    if abs(-power/record["power"]-1.0) > tolerance:
        sys.exit(f"{name}: spectrum integrates to {-power:.4e}, header power is {record['power']:.4e}")
    umin, umax = record["umin"], record["umax"]
    intensityMean = umin if umax == umin else np.log(umax/umin)/(1.0/umin-1.0/umax)
    if abs(intensityMean/record["intensityMean"]-1.0) > tolerance:
        sys.exit(f"{name}: <U> in header is {record['intensityMean']:.4e}, expected {intensityMean:.4e}")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: draineLi2007ToHDF5.py DL07spec.tgz")
    pattern = re.compile(r"^U([0-9.e]+)/U\1_([0-9.e]+)_([A-Za-z0-9._]+)\.txt$")
    records = {}
    with tarfile.open(sys.argv[1], "r:gz") as tar:
        for member in tar.getmembers():
            match = pattern.match(member.name)
            if not member.isfile() or not match:
                continue
            record = parse(tar.extractfile(member).read().decode(), member.name)
            if record is not None:
                records[(match.group(3), float(match.group(1)), float(match.group(2)))] = record

    # The umin grid: values for which every model has a single-intensity spectrum and every power-law spectrum.
    intensityMinimum = sorted({key[1] for key in records})
    intensityMinimum = np.array([umin for umin in intensityMinimum
                                 if all((model[0], umin, umax) in records for model in models for umax in (umin, *intensityMaximum))])
    wavelength       = records[(models[0][0], intensityMinimum[0], intensityMinimum[0])]["wavelength"]

    path = os.path.join(os.environ.get("GALACTICUS_DATA_PATH", "."), "static", "dust", "emission")
    os.makedirs(path, exist_ok=True)
    fileHDF5 = os.path.join(path, "draineLi2007.hdf5")
    with h5py.File(fileHDF5, "w") as file:
        file.create_dataset("wavelength"      , data=wavelength      )
        file.create_dataset("intensityMinimum", data=intensityMinimum)
        file.create_dataset("intensityMaximum", data=intensityMaximum)
        file["wavelength"      ].attrs["units"      ] = "Angstroms"
        file["intensityMinimum"].attrs["description"] = "minimum starlight intensity, U, in units of the local interstellar radiation field of Mathis, Mezger & Panagia (1983)"
        file["intensityMaximum"].attrs["description"] = "maximum starlight intensity of the power-law distributions"
        for nameFile, label, fractionMassPAH, massDustPerHydrogen in models:
            single       = np.zeros((                      intensityMinimum.size, wavelength.size))
            powerLaw     = np.zeros((intensityMaximum.size, intensityMinimum.size, wavelength.size))
            powerSingle  = np.zeros((                      intensityMinimum.size                 ))
            powerPowered = np.zeros((intensityMaximum.size, intensityMinimum.size                 ))
            for i, umin in enumerate(intensityMinimum):
                for j, umax in enumerate((umin, *intensityMaximum)):
                    name   = f"U{umin}_{umax}_{nameFile}"
                    record = records[(nameFile, umin, umax)]
                    if not np.array_equal(record["wavelength"], wavelength):
                        sys.exit(f"{name}: wavelength grid differs")
                    check(record, name)
                    if j == 0:
                        single     [   i] = record["emissivity"]
                        powerSingle[   i] = record["power"     ]
                    else:
                        powerLaw    [j-1, i] = record["emissivity"]
                        powerPowered[j-1, i] = record["power"     ]
            group = file.create_group(label)
            group.create_dataset("emissivitySingle"  , data=single      , compression="gzip", compression_opts=9)
            group.create_dataset("emissivityPowerLaw", data=powerLaw    , compression="gzip", compression_opts=9)
            group.create_dataset("powerSingle"       , data=powerSingle )
            group.create_dataset("powerPowerLaw"     , data=powerPowered)
            for dataset in ("emissivitySingle", "emissivityPowerLaw"):
                group[dataset].attrs["description"] = "nu dP/dnu, the power radiated per unit logarithmic frequency per H nucleon"
                group[dataset].attrs["units"      ] = "erg/s"
            for dataset in ("powerSingle", "powerPowerLaw"):
                group[dataset].attrs["description"] = "power radiated per H nucleon"
                group[dataset].attrs["units"      ] = "erg/s"
            group.attrs["name"               ] = nameFile
            group.attrs["fractionMassPAH"    ] = fractionMassPAH
            group.attrs["massDustPerHydrogen"] = massDustPerHydrogen
        file.attrs["source"] = "Draine & Li (2007), ApJ, 657, 810; spectra from irem4/DL07spec.tgz, and Table 3"
        file.attrs["url"   ] = "https://www.astro.princeton.edu/~draine/dust/irem.html"
    print(f"draineLi2007: wavelength {wavelength.size} ({wavelength.min():.4g}-{wavelength.max():.4g} A), "
          f"intensityMinimum {intensityMinimum.size} ({intensityMinimum.min()}-{intensityMinimum.max()}), "
          f"intensityMaximum {intensityMaximum.size}, models {len(models)}")


if __name__ == "__main__":
    main()
