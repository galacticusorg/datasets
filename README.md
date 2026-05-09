# Galacticus datasets

Static and dynamically-generated data files used by [Galacticus](https://github.com/galacticusorg/galacticus).

## Layout

| Path        | Purpose |
| ----------- | ------- |
| `static/`   | Curated input data shipped with the repo. Read-only at run time. |
| `dynamic/`  | Populated at run time by Galacticus (and by analysis/post-processing) with derived datasets. Git-ignored. |
| `scripts/`  | Maintenance scripts run from CI (e.g. workflow status reporter). |
| `.github/`  | CI workflows that update and validate the Local Group satellite database, diff HDF5 files in PRs, and report workflow status to Slack. |

`static/` is organised by physics topic — `darkMatter/`, `cooling/`, `chemicalState/`, `stellarPopulations/`, `filters/`, `observations/`, `surveyGeometry/`, etc. Each subdirectory holds the data files (typically HDF5 or XML) referenced from the corresponding Galacticus modules.

## How Galacticus consumes it

Galacticus locates this tree via the `GALACTICUS_DATA_PATH` environment variable (or its build-time default). At run time it reads files from `static/` and writes any derived/cached data into `dynamic/`. Renaming files or directories here therefore requires a matching change in the [`galacticus`](https://github.com/galacticusorg/galacticus) repository.

## Working with the data

Two git diff drivers are configured in `.gitattributes` so that PRs show meaningful diffs of binary HDF5 and large XML files:

* `*.hdf5` → `h5diff -r -c` (requires `hdf5-tools`)
* `*.xml`  → `xdiff` (Galacticus' XML differ)

Enable them locally with:

```sh
git config --local include.path ../.gitconfig
```

The `HDF5-Diff` GitHub Actions workflow runs the same diff on PRs and posts the output as a comment, with links to a hosted HDF5 viewer.

## Local Group satellite database

`static/observations/localGroup/localGroupSatellites.xml` is validated against `localGroupSatellites.xsd` on every PR (`Validate-Local-Group-DB` workflow), and its bibliographic references are refreshed weekly against NASA ADS by `Update-Local-Group-DB`. Python dependencies for both workflows are pinned in `static/observations/localGroup/requirements.txt`.

## License and upstream sources

This repository is released under the MIT License (see `LICENSE`).

A substantial fraction of the files under `static/` are *redistributions* of data originally published elsewhere — for example BC2003 SSP spectra, ULTRAVISTA completeness tables (see `static/surveyGeometry/ULTRAVISTA/README`), the Wilms (2000) ISM absorption table, McConnell &amp; Ma black-hole kinematics, FSPS source patches, etc. The MIT License grant in this repository covers the curated *compilation* and any genuinely-original content (such as `localGroupSatellites.xml` and the Galacticus-generated cosmological parameter files); it does **not** override whatever terms each upstream source attached to the underlying data. Users who plan to redistribute individual files should consult the original publication or data release for terms applying to that file.
