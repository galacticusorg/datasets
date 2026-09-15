#!/usr/bin/env python3
"""Reduce the PAH emission basis spectra of Richie & Hensley (2026) to a compact table for Galacticus.

In the single-photon approximation the emission of a PAH population is linear in the radiation heating it: each photon
absorbed produces the emission spectrum of a single absorption, independently of any other. Richie & Hensley (2026;
arXiv:2510.16861) publish those "basis" spectra, p(a, lambda_abs, lambda_em), for neutral and ionized PAHs of 59 sizes
(3.5-100 Angstroms), absorbing photons in 474 bins of fractional width 1% from 912 Angstroms to 10.1 microns, emitting
at 8944 wavelengths from 0.1 micron to 1 cm. They are distributed with the code pah_spec
(https://github.com/helenarichie/pah_spec, MIT license) as two 1 GB files, basis_ion.h5 and basis_neu.h5, at
https://doi.org/10.7910/DVN/LUUXEJ (CC0).

A galaxy model needs only the emission of a whole PAH population per unit energy absorbed in each bin. That is
computed here, once, for each of the size distributions (small, standard, large) and ionization functions (low,
standard, high) of Draine et al. (2021), as implemented by pah_spec:

  * The energy a population absorbs in bin i is shared among sizes and charge states in proportion to
    dn(a) C_abs(a, lambda), with C_abs averaged over the bin (the radiation field is taken as flat across a 1% bin).
  * Each basis spectrum is normalized to unit power, weighted by that share, and summed over sizes and charge states.
  * The result is resampled, conserving energy, onto emission bins of resolution R=500 from 1 to 20 microns and R=100
    elsewhere, and stored as the bin average of lambda p_lambda per unit ln(lambda), which sums (times the bin widths
    in ln lambda) to one for every absorbed-photon bin.

Energy must also be shared between PAHs and the other grains which absorb the same light. So the PAH absorption
opacity per H of each population, sum_a dn(a) C_abs(a, lambda), is stored for each bin, along with the absorption
opacity per H of astrodust from the Astrodust+PAH model of Hensley & Draine (2023; https://doi.org/10.7910/DVN/3B6E6S,
CC0). The astrodust tabulation begins at 0.1 micron, so below that (912-1000 Angstroms) it is extrapolated as a power
law in wavelength fitted to its first two points. The PAH opacity is taken from pah_spec rather than from the "PAH"
component of Astrodust+PAH, so that the PAHs which absorb are the ones whose emission is tabulated; the two differ,
pah_spec's population absorbing about two thirds as much in the ultraviolet.

Usage:

    richieHensley2026ToHDF5.py BASIS_DIRECTORY PAH_SPEC_DATA_DIRECTORY ASTRODUST_FITS [--no-validate]

where BASIS_DIRECTORY holds basis_ion.h5 and basis_neu.h5, PAH_SPEC_DATA_DIRECTORY is the "data" directory of a pah_spec
checkout (for its absorption cross sections), and ASTRODUST_FITS is astrodust+PAH_MW_RV3.1.fits. pah_spec must be
importable. The layout written, to $GALACTICUS_DATA_PATH/static/dust/emission/richieHensley2026.hdf5, is:

    /wavelengthAbsorbedEdges     (475)                  Angstroms; bin i spans edges i to i+1
    /wavelengthEmittedEdges      (nEmitted+1)           Angstroms
    /opacityAbsorptionAstrodust  (474)                  cm^2 per H, averaged over each absorbed bin
    /<population>/opacityAbsorptionPAH  (474)           cm^2 per H, averaged over each absorbed bin
    /<population>/emission       (474, nEmitted)        mean lambda p_lambda per ln(lambda) per unit absorbed energy

with populations named size{Small,Standard,Large}Ionization{Low,Standard,High}.

Unless --no-validate is given, the reduction is checked against pah_spec itself: the spectrum of the standard
population heated by the U=1 interstellar radiation field of Mathis, Mezger & Panagia (1983), computed by pah_spec from
the full basis spectra, is compared with the same spectrum rebuilt from this table, bin by bin.

Andrew Benson, with assistance from Claude.
"""
import os
import sys

import astropy.units as u
import h5py
import numpy as np
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid, trapezoid

import pah_spec

micronsToAngstroms = 1.0e4
speedOfLight       = 2.99792e10  # cm/s, as used by pah_spec
widthAbsorbed      = 0.01        # fractional width of the absorbed-photon bins
samplesPerBin      = 8           # samples per absorbed bin when averaging cross sections
sizes              = {"Small": "sm", "Standard": "st", "Large": "lg"}
ionizations        = {"Low": "lo", "Standard": "st", "High": "hi"}


def emissionEdges():
    """Emission bin edges in microns: R=100 from 0.1 to 1 micron, R=500 from 1 to 20 microns, R=100 from 20 microns to 1 cm."""
    edges = []
    for lower, upper, resolution in ((0.1, 1.0, 100.0), (1.0, 20.0, 500.0), (20.0, 1.0e4, 100.0)):
        count = int(np.ceil(np.log(upper/lower)*resolution))
        segment = np.exp(np.linspace(np.log(lower), np.log(upper), count+1))
        edges.append(segment if not edges else segment[1:])
    return np.concatenate(edges)


def main():
    arguments = [argument for argument in sys.argv[1:] if not argument.startswith("--")]
    validate  = "--no-validate" not in sys.argv[1:]
    if len(arguments) != 3:
        sys.exit("usage: richieHensley2026ToHDF5.py BASIS_DIRECTORY PAH_SPEC_DATA_DIRECTORY ASTRODUST_FITS [--no-validate]")
    directoryBasis, directoryData, fileAstrodust = arguments

    # Read the grids and the basis spectra, normalizing each to unit power as it is read.
    basis = {}
    for charge in ("ion", "neu"):
        with h5py.File(os.path.join(directoryBasis, f"basis_{charge}.h5"), "r") as file:
            wavelengthAbsorbed = file["lambda_abs" ][:]
            wavelengthEmitted  = file["lambda_em"  ][:]
            radii              = file["grain_sizes"][:]
            spectra            = file["basis_spectra"][:]
        power = trapezoid(spectra, x=wavelengthEmitted, axis=2)
        if np.any(power <= 0.0):
            sys.exit(f"basis_{charge}.h5: a basis spectrum has no power")
        spectra /= power[:, :, np.newaxis].astype(spectra.dtype)
        basis[charge] = spectra
        print(f"basis_{charge}.h5: {spectra.shape} read and normalized", flush=True)
    if not np.allclose(wavelengthAbsorbed[1:]/wavelengthAbsorbed[:-1], 1.0+widthAbsorbed, rtol=1.0e-9, atol=0.0):
        sys.exit("absorbed-photon bins are not contiguous with the expected fractional width")
    edgesAbsorbed = np.append(wavelengthAbsorbed, wavelengthAbsorbed[-1]*(1.0+widthAbsorbed))

    # Absorption cross sections from pah_spec, averaged over each absorbed bin. The PahSpec object is built on the small
    # test basis shipped with pah_spec, which shares the grain sizes, to avoid holding the full basis twice.
    helper = pah_spec.PahSpec(basis_dir=os.path.join(directoryData, "test_data", "basis_spectra_low_res"), internal_data_dir=directoryData)
    if not np.allclose(helper.grain_sizes.to(u.AA).value, radii, rtol=1.0e-9, atol=0.0):
        sys.exit("grain sizes of the full basis differ from those of pah_spec's test basis")
    fractions  = (np.arange(samplesPerBin)+0.5)/samplesPerBin
    samples    = (wavelengthAbsorbed[:, np.newaxis]*(1.0+widthAbsorbed*fractions[np.newaxis, :])).ravel()
    crossIon, crossNeutral = helper.calc_c_abs(samples*u.um, helper.grain_sizes)
    crossSection = {
        "ion": np.asarray(crossIon    .to(u.cm**2).value).reshape(radii.size, wavelengthAbsorbed.size, samplesPerBin).mean(axis=2),
        "neu": np.asarray(crossNeutral.to(u.cm**2).value).reshape(radii.size, wavelengthAbsorbed.size, samplesPerBin).mean(axis=2),
    }

    # Astrodust absorption opacity per H, averaged over each absorbed bin, extrapolated below the tabulation.
    with fits.open(fileAstrodust) as file:
        extinction = file["EXTINCTION"].data
        scattering = file["SCATTERING"].data
    wavelengthAstrodust = extinction[:, 0]
    opacityAstrodust    = extinction[:, 1]-scattering[:, 1]
    slope               = np.log(opacityAstrodust[1]/opacityAstrodust[0])/np.log(wavelengthAstrodust[1]/wavelengthAstrodust[0])
    sampled             = np.exp(np.interp(np.log(samples), np.log(wavelengthAstrodust), np.log(opacityAstrodust)))
    below               = samples < wavelengthAstrodust[0]
    sampled[below]      = opacityAstrodust[0]*(samples[below]/wavelengthAstrodust[0])**slope
    opacityAstrodustBins = sampled.reshape(wavelengthAbsorbed.size, samplesPerBin).mean(axis=1)

    # Emission bins, and the cumulative integrals used to resample onto them.
    edgesEmitted = emissionEdges()
    widthsLog    = np.log(edgesEmitted[1:]/edgesEmitted[:-1])

    def resample(spectrum):
        """Resample spectra normalized to unit power (last axis on the native grid) to mean lambda p_lambda per ln lambda."""
        cumulative = cumulative_trapezoid(spectrum, x=wavelengthEmitted, axis=-1, initial=0.0)
        atEdges    = np.stack([np.interp(edgesEmitted, wavelengthEmitted, row) for row in cumulative.reshape(-1, wavelengthEmitted.size)])
        return (np.diff(atEdges, axis=-1)/widthsLog).reshape(spectrum.shape[:-1]+(widthsLog.size,))

    populations = {}
    path = os.path.join(os.environ.get("GALACTICUS_DATA_PATH", "."), "static", "dust", "emission")
    os.makedirs(path, exist_ok=True)
    fileHDF5 = os.path.join(path, "richieHensley2026.hdf5")
    with h5py.File(fileHDF5, "w") as file:
        file.create_dataset("wavelengthAbsorbedEdges"   , data=edgesAbsorbed*micronsToAngstroms)
        file.create_dataset("wavelengthEmittedEdges"    , data=edgesEmitted *micronsToAngstroms)
        file.create_dataset("opacityAbsorptionAstrodust", data=opacityAstrodustBins)
        file["wavelengthAbsorbedEdges"   ].attrs["units"      ] = "Angstroms"
        file["wavelengthEmittedEdges"    ].attrs["units"      ] = "Angstroms"
        file["opacityAbsorptionAstrodust"].attrs["units"      ] = "cm^2 per H"
        file["opacityAbsorptionAstrodust"].attrs["description"] = "absorption opacity of astrodust (Hensley & Draine 2023), averaged over each absorbed-photon bin; extrapolated as a power law below 0.1 micron"
        for labelSize, size in sizes.items():
            for labelIonization, ionization in ionizations.items():
                label = f"size{labelSize}Ionization{labelIonization}"
                numberNeutral, numberIon = pah_spec.calc_dn(helper.grain_sizes, size_dist=size, ion_frac=ionization)
                numberNeutral = np.asarray(getattr(numberNeutral, "value", numberNeutral))
                numberIon     = np.asarray(getattr(numberIon    , "value", numberIon    ))
                weightNeutral = numberNeutral[:, np.newaxis]*crossSection["neu"]
                weightIon     = numberIon    [:, np.newaxis]*crossSection["ion"]
                opacityPAH    = (weightNeutral+weightIon).sum(axis=0)
                combined      = (np.einsum("ai,aie->ie", (weightNeutral/opacityPAH).astype(np.float32), basis["neu"])
                                 +np.einsum("ai,aie->ie", (weightIon    /opacityPAH).astype(np.float32), basis["ion"]))
                emission      = resample(combined.astype(np.float64))
                sums          = (emission*widthsLog).sum(axis=1)
                if np.max(np.abs(sums-1.0)) > 1.0e-3:
                    sys.exit(f"{label}: resampled emission does not conserve energy (worst {np.max(np.abs(sums-1.0)):.2e})")
                group = file.create_group(label)
                group.create_dataset("opacityAbsorptionPAH", data=opacityPAH)
                group.create_dataset("emission"            , data=emission.astype(np.float32), compression="gzip", compression_opts=9)
                group["opacityAbsorptionPAH"].attrs["units"      ] = "cm^2 per H"
                group["emission"            ].attrs["description"] = "mean lambda p_lambda per unit ln(lambda) in each emitted bin, per unit energy absorbed in each absorbed-photon bin"
                group.attrs["sizeDistribution"   ] = size
                group.attrs["ionizationFunction"] = ionization
                populations[label] = (opacityPAH, emission)
                print(f"{label}: energy conserved to {np.max(np.abs(sums-1.0)):.1e}; PAH/astrodust absorption at 0.1 micron "
                      f"{np.interp(0.1, wavelengthAbsorbed, opacityPAH/opacityAstrodustBins):.3f}", flush=True)
        file.attrs["source"] = "Richie & Hensley (2026), arXiv:2510.16861; basis spectra doi:10.7910/DVN/LUUXEJ (pah_spec v1.0.0); astrodust absorption from Hensley & Draine (2023), doi:10.7910/DVN/3B6E6S"
        file.attrs["url"   ] = "https://github.com/helenarichie/pah_spec"
    print(f"wrote {fileHDF5}: {wavelengthAbsorbed.size} absorbed bins, {widthsLog.size} emitted bins, {len(populations)} populations", flush=True)

    if not validate:
        return
    # Validation: pah_spec's own spectrum of the standard population in the U=1 MMP83 field, from the full basis.
    del basis
    reference = pah_spec.PahSpec(basis_dir=directoryBasis, internal_data_dir=directoryData)
    neutral, ion = reference.generate_spectrum()
    spectrumReference = (neutral+ion).to(u.erg/u.s/u.cm).value
    wavelengthsField, energyDensityField = reference.wavelength_u_arr.to(u.um).value, reference.u_lambda_arr.to(u.erg/u.cm**3/u.um).value
    del reference
    powerReference = trapezoid(spectrumReference, x=wavelengthEmitted*1.0e-4)
    # Rebuild from the table: the energy absorbed per H in each bin, times the normalized emission.
    opacityPAH, emission = populations["sizeStandardIonizationStandard"]
    fieldSampled = np.interp(samples, wavelengthsField, energyDensityField).reshape(wavelengthAbsorbed.size, samplesPerBin).mean(axis=1)
    absorbed     = speedOfLight*fieldSampled*opacityPAH*np.diff(edgesAbsorbed)   # erg/s per H
    rebuilt      = absorbed @ emission                                           # lambda p_lambda per ln lambda, erg/s per H
    # Resample pah_spec's spectrum onto the same bins, as power per ln lambda.
    cumulative   = cumulative_trapezoid(spectrumReference*1.0e-4, x=wavelengthEmitted, initial=0.0)
    referenceBins = np.diff(np.interp(edgesEmitted, wavelengthEmitted, cumulative))/widthsLog
    powerRebuilt = (rebuilt*widthsLog).sum()
    band         = (edgesEmitted[:-1] >= 3.0) & (edgesEmitted[1:] <= 20.0) & (referenceBins > 1.0e-3*referenceBins.max())
    worst        = np.max(np.abs(rebuilt[band]/referenceBins[band]-1.0))
    print(f"validation (standard population, MMP83 U=1): total power rebuilt/pah_spec = {powerRebuilt/powerReference:.5f}; "
          f"worst bin at 3-20 microns differs by {worst:.2e}", flush=True)
    if abs(powerRebuilt/powerReference-1.0) > 5.0e-3 or worst > 2.0e-2:
        sys.exit("validation failed: the reduced table does not reproduce pah_spec")


if __name__ == "__main__":
    main()
