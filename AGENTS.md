# AGENTS.md

Guidance for AI coding agents (Claude Code, Codex, Copilot, Gemini CLI, Cursor,
…) working in this repository; humans may find it useful too. Claude Code loads
it through `.claude/CLAUDE.md`, which imports it. `README.md` describes the
layout, the diff drivers and the Local Group satellite database; this file
adds what an agent needs before changing anything.

This repository holds the data files that
[Galacticus](https://github.com/galacticusorg/galacticus) reads at run time,
through `GALACTICUS_DATA_PATH`. That repository's own
[`AGENTS.md`](https://github.com/galacticusorg/galacticus/blob/master/AGENTS.md)
applies whenever you build or run Galacticus to check a data change.

## Cloning

The repository is large, and most tasks touch one directory. Clone without
blobs and check out only what you need:

```sh
git clone --filter=blob:none --sparse https://github.com/galacticusorg/datasets.git
git sparse-checkout add static/<topic>
```

In such a clone, avoid commands that read every blob (`git ls-tree -l`,
`git log -p` or `git grep` over the whole tree): each one fetches the missing
blobs, one round trip at a time.

## Changing data

- **`static/` is read by Galacticus code by name.** Renaming or moving a file or
  directory needs a matching change in the `galacticus` repository; open both
  together and link them. `dynamic/` is written at run time and git-ignored:
  never commit anything in it.
- **A changed file must still be read correctly.** Run the Galacticus code (a
  test in `testSuite/`, or a small model) that reads the file, against this
  checkout as `GALACTICUS_DATA_PATH`, and say in the PR what you ran.
- **Keep generated files reproducible.** A file produced by a script records
  its provenance in a comment block at the top (see
  `static/stellarAstrophysics/Supernovae_Type_Ia_Yields_Gronow2021_*.xml`): the
  science source and URL, any transcription source, and the script that
  generated it. Such scripts live in
  [`galacticusDevTools`](https://github.com/galacticusorg/galacticusDevTools)
  (`dataPreparation/`); change the script and regenerate, rather than editing
  the output by hand.
- **`static/observations/localGroup/localGroupSatellites.xml`** is validated
  against its XSD on every PR, and its references are refreshed weekly by a
  workflow. Validate locally before pushing, with the dependencies in that
  directory's `requirements.txt`.
- Keep XML files ASCII (write `&amp;`, `&lt;` and so on), and use US English
  spelling in descriptions and comments.

## Credit for published data

Reusing a published academic table needs **credit, not a license**. Record it
in the file:

- one file-level `<source>` element naming the specific table or figure (for
  example `Heger &amp; Woosley (2002; ApJ; 567; 532; Figure 1)`), and
- one `<url>` element, preferring a direct link to the data (a journal's
  machine-readable table, VizieR, Zenodo) over an ADS abstract page.

`static/stellarAstrophysics/` shows the convention. For HDF5 files, put the
same information in file-level attributes, following a sibling file in the same
directory.

Raise licensing only when the source attaches real conditions (a CC-BY license
on Zenodo, say, or terms on a data release page): then record those terms, and
ask the maintainer before adding the file. The MIT license of this repository
covers the compilation and original content only, not the terms of upstream
data (see *License and upstream sources* in `README.md`).

## Commits and pull requests

- Use [Conventional Commits](https://www.conventionalcommits.org)
  (`type(scope): summary`, with the types Galacticus allows: `fix`, `feat`,
  `build`, `chore`, `ci`, `docs`, `style`, `test`, `refactor`, `perf`,
  `revert`, `clean`).
- Work on a branch and open a pull request; don't push to `master`. A
  maintainer reviews and merges. The `HDF5-Diff` workflow posts a diff of
  changed HDF5 files on the PR: check that it shows only the change you meant.
- Describe in the commit message what an AI tool contributed and how the result
  was verified, as Galacticus' `CONTRIBUTING.md` asks. All AI-generated work
  must be verified by a human before submission.
