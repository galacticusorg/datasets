#!/usr/bin/env python3
"""Convert the Dale & Helou (2002) infrared spectral energy distributions of star-forming galaxies to HDF5.

The templates are distributed at http://physics.uwyo.edu/~ddale/research/seds/seds.html as two ASCII files:

    spectra.dat   column 1 the wavelength in microns; columns 2-65 log10(nu f_nu) for each of 64 values of alpha
    alpha.dat     column 1 those 64 values of alpha; column 2 log10 of the IRAS f_nu(60um)/f_nu(100um) ratio

Here alpha is the exponent of the distribution of dust mass over the intensity U of the heating radiation field,
dM(U) ~ U^-alpha dU for 0.3 <= U <= 1e5 (Dale et al. 2001, ApJ, 549, 215). The nominal units of nu f_nu are W per H
atom, but the author states that this absolute scale has never been used; the templates are shapes, and the model
normalizes them.

Usage:

    daleHelou2002ToHDF5.py spectra.dat alpha.dat

The layout written, to $GALACTICUS_DATA_PATH/static/dust/emission/daleHelou2002.hdf5, is:

    /wavelength   (nWavelength)          Angstroms, ascending
    /alpha        (nAlpha)               exponent of the dust mass distribution over heating intensity
    /luminosity   (nAlpha, nWavelength)  nu L_nu, arbitrary units (nominally W per H atom)

The tabulation is not all dust emission. Checked on conversion, and recorded as attributes so that the model can not
lose it:

  * Below 3 microns the shape of every template is the same to 0.002 dex whatever alpha, so it is a fixed component
    independent of the dust heating (presumably stellar), carrying 5% (alpha=0.0625) to 38% (alpha=4) of the tabulated
    power. The paper presents its spectra from 3 microns.
  * Beyond about 1 mm the templates are extended to radio wavelengths using the far-infrared--radio correlation
    (synchrotron and free-free emission), which is not emission from dust.

So `wavelengthDustMinimum` and `wavelengthDustMaximum` (3 and 1100 microns, the range the author quotes for these
spectra) bound the dust emission. The conversion itself is faithful: no points are added, removed, or altered.

Andrew Benson, with assistance from Claude.
"""
import os
import sys

import h5py
import numpy as np

micronsToAngstroms    = 1.0e4
wavelengthDustMinimum = 3.0e0*micronsToAngstroms
wavelengthDustMaximum = 1.1e3*micronsToAngstroms


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: daleHelou2002ToHDF5.py spectra.dat alpha.dat")
    fileSpectra, fileAlpha = sys.argv[1:3]

    spectra = np.loadtxt(fileSpectra)
    alpha   = np.loadtxt(fileAlpha  )[:, 0]
    if spectra.ndim != 2 or spectra.shape[1] != alpha.size+1:
        sys.exit(f"{fileSpectra}: expected {alpha.size+1} columns, found shape {spectra.shape}")
    if not np.all(np.diff(alpha) > 0.0):
        sys.exit(f"{fileAlpha}: alpha is not strictly increasing")
    wavelength = spectra[:, 0]*micronsToAngstroms
    if not np.all(np.diff(wavelength) > 0.0):
        sys.exit(f"{fileSpectra}: wavelength is not strictly increasing")
    luminosity = 10.0**spectra[:, 1:].T

    # Confirm that the component below the dust range is independent of alpha, as described above; if a future release
    # changes that, the dust range recorded here may no longer be right.
    blue  = wavelength < wavelengthDustMinimum
    shape = np.log10(luminosity[:, blue])-np.log10(luminosity[:, blue][:, -1:])
    if np.ptp(shape, axis=0).max() > 0.01:
        sys.exit(f"{fileSpectra}: the component below 3 microns now depends on alpha; revisit the dust range")

    path = os.path.join(os.environ.get("GALACTICUS_DATA_PATH", "."), "static", "dust", "emission")
    os.makedirs(path, exist_ok=True)
    fileHDF5 = os.path.join(path, "daleHelou2002.hdf5")
    with h5py.File(fileHDF5, "w") as file:
        file.create_dataset("wavelength", data=wavelength)
        file.create_dataset("alpha"     , data=alpha     )
        file.create_dataset("luminosity", data=luminosity)
        file["wavelength"].attrs["units"      ] = "Angstroms"
        file["alpha"     ].attrs["description"] = "exponent of the distribution of dust mass over heating intensity, dM ~ U^-alpha dU"
        file["luminosity"].attrs["description"] = "nu L_nu, in arbitrary units (nominally W per H atom)"
        file.attrs["source"               ] = "Dale & Helou (2002), ApJ, 576, 159; spectra.dat and alpha.dat"
        file.attrs["url"                  ] = "http://physics.uwyo.edu/~ddale/research/seds/seds.html"
        file.attrs["wavelengthDustMinimum"] = wavelengthDustMinimum
        file.attrs["wavelengthDustMaximum"] = wavelengthDustMaximum
        file.attrs["note"                 ] = ("Only wavelengthDustMinimum <= wavelength <= wavelengthDustMaximum is dust emission: shortward is an "
                                               "alpha-independent (presumably stellar) component, longward a radio extension from the "
                                               "far-infrared--radio correlation.")
    print(f"daleHelou2002: wavelength {wavelength.size} ({wavelength.min():.4g}-{wavelength.max():.4g} A), alpha {alpha.size}, "
          f"dust range points {np.count_nonzero((wavelength >= wavelengthDustMinimum) & (wavelength <= wavelengthDustMaximum))}")


if __name__ == "__main__":
    main()
