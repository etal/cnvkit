# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CNVkit is a command-line toolkit and Python library for detecting copy number variants and alterations genome-wide from high-throughput sequencing data. It provides both a CLI interface and Python API for genomic analysis workflows.

**Supported Python versions:** 3.11+ (tested on 3.11-3.14)
**Minimum versions aligned with Ubuntu 26.04 LTS (Resolute)**

## Development Workflow

- **Bug fixes and new features**: Write a failing test first, then implement.
- **Edge cases**: Before finishing, verify behavior for empty inputs, NaN/missing values, and single-element arrays.
- **NaN weight safety**: `np.average(values, weights=w)` and `numpy.sum(w)` propagate NaN; use `np.nansum` for weight sums and filter `~np.isnan(wt)` before `np.average`. Note that `pandas.Series.sum()` skips NaN by default — prefer explicit `np.nansum` for clarity.
- **User-facing changes**: Update the relevant docs in `doc/*.rst`.
- **Clinical impact**: When reviewing changes, consider whether the changeset alters numerical output or output file formats (`.cnr`, `.cns`, `.cnn`, SEG, VCF). Flag any such changes explicitly, as downstream clinical pipelines may depend on exact output stability.

## Development Commands

See `doc/development.rst` for the full developer guide (environment
setup, pre-commit hooks, code style, testing matrix, PR process, Docker
release flow). The notes below are the AI-session-specific quirks not
already in that guide.

**Conda env quirk:** `conda run -n cnvkit <command>` is unreliable in
this repo (the dev tools sometimes fail to find their dependencies).
Use `conda activate cnvkit && <command>` instead — including inside
scripts and parallel agent runs.

**`# type: ignore` requirement:** `pyproject.toml` sets
`enable_error_code = ["ignore-without-code"]`, so every `# type: ignore`
must specify the error code (e.g. `# type: ignore[return-value]`). Bare
`# type: ignore` fails mypy.

**Common type-ignore patterns in this codebase** (the footguns mypy
flags most often here):
- `tabio.read()` returns union types — `# type: ignore[return-value]`
  at call sites.
- Pandas operations often return `Any` — `# type: ignore[no-any-return]`.
- Closures don't narrow types in mypy — `# type: ignore[index]` or a
  local `assert`.
- Parameters typed `param: None = None` produce unreachable blocks —
  use `Optional[Type] = None` (or `Type | None = None`) instead.
- Generator functions must use `-> Generator[YieldType, SendType,
  ReturnType]`, not `-> Iterator` or `-> None`.
- When a variable changes type (e.g. `str` → `list[int]`), rename it to
  avoid shadowing (e.g. `copies` → `copy_strs`).
- `numpy.bool_` is not assignable to `bool` in mypy — use
  `# type: ignore[assignment]` or widen parameters to
  `bool | bool_ | None`.

## Architecture

### Core Structure
- **`cnvlib/`** - Main Python package
  - `commands.py` - CLI definitions and API functions (`_cmd_*` for arg parsing, `do_*` for logic)
  - `cnvkit.py` - CLI entry point that routes to commands
  - `core.py` - Core data structures and utilities
  - `segmentation/` - Segmentation algorithms (CBS, HMM, etc.)
  - `batch.py`, `segment.py`, etc. - Individual command implementations
  - `cmdutil.py`, `params.py` - Utility functions
  - `plots.py`, `diagram.py`, `heatmap.py`, `scatter.py` - Visualization
  - `importers.py`, `export.py` - Data import/export
  - `cnary.py` - CopyNumArray (extends GenomicArray with log2 ratios and gene names)
  - `vary.py` - VariantArray (extends GenomicArray with variant allele data)

- **`skgenome/`** - Genomic data handling library (part of CNVkit but decoupled)
  - `gary.py` - GenomicArray class for genomic interval data (wraps pandas DataFrame)
  - `tabio/` - File I/O for BED, GFF, VCF, SEG, Picard, CNVkit formats (.cnn/.cnr/.cns), bedGraph
  - `chromnames.py` - Chromosome-name classification (autosome/sex/mito/alt-contig detection, arabic + Roman numerals, X/Y label inference)
  - `genomebuild.py` - Reference assembly metadata (PAR coordinates) as `GenomeBuild` value objects, with `get_genome_build()` lookup
  - `chromsort.py` - Chromosome-name sort keys

GenomicArray uses `typing.Self` (PEP 673) so that methods like `.copy()`, `.concat()`, and `.as_dataframe()` preserve the subclass type through type checking.

### Chromosome-name handling

**All chromosome-name classification goes through `skgenome.chromnames`** —
do NOT add inline regexes or `.startswith("chr")` / `chromosome.iat[0]`
heuristics. The classifier is *context-aware* (it inspects the whole
chromosome set, not one name), because `chrX` is a sex chromosome in
human but autosome 10 in yeast. See `doc/sex.rst` "Non-human and
Roman-numeral genomes" for the user-facing behavior.

API points worth knowing:
- `CopyNumArray.chr_x_label` / `chr_y_label` return `str | None`. `None`
  means "no sex chromosome detected in this assembly" — callers must
  handle `None` (the `chr_*_filter` methods return all-False in that
  case rather than crashing).
- `GenomicArray.autosomes()` falls back to returning the whole array
  (with a warning) when no autosomes are recognized. Be permissive on
  unfamiliar assemblies; don't silently drop data.
- PAR coordinates live in `skgenome.genomebuild`;
  `cnvlib.params.PSEUDO_AUTSOMAL_REGIONS` is a back-compat re-export.

### Sex inference

`doc/sex.rst` is the source of truth for the math (ratio-of-residuals
maleness ratios + AND-gate, VCF heterozygous-SNP confirmer, target /
antitarget reconciliation in `do_reference`). The rules below are
project invariants that, if broken, silently regress the design
without breaking tests:

- **`verify_sample_sex` is the canonical resolver.** Don't reinvent
  `is_female_default(guess_xx(...))` inline; route through it so user
  `--sample-sex` and the VCF het-density confirmer both apply.
- **Honest `None` propagates on inference paths**; concrete `bool` only
  at decision consumers via `is_female_default`. `do_reference`'s
  target/antitarget reconciliation needs the honest `None` to work.
- **Reporting commands stay honest**: `do_sex` reports `Unknown` for an
  undeterminable sample; don't collapse `None` → `Female` at the report
  layer.
- **VCF het confirmer is one-way (male → female only).** Absent chrX
  hets is non-evidence — true haploid X and "too few SNPs" are
  indistinguishable.
- Do NOT re-introduce a chi-square or multiplicative `combined_score`
  on top of the median; the median is already the robust quantity, and
  wrapping it in a chi-square brought sample-size dependence that
  produced `.cnr` vs `.cns` inconsistency (#785).

### File Formats
- `.cnn` - Coverage/reference data
- `.cnr` - Copy number ratio data
- `.cns` - Segmented copy number data

## Dependencies

Core dependencies are declared in `requirements/core.txt`; `min.txt` pins exact minimums for compatibility testing.

## Packaging

Packaging is deliberately conservative. CNVkit ships across PyPI, conda, Docker, and Galaxy, and that matrix is far more fragile than the numerical core, so the project avoids formal `package-data`:

- Small, code-adjacent data stays inline (e.g. PAR coordinates as Python literals in `skgenome.genomebuild`; the CBS R script as a string in `cnvlib/segmentation/cbs.py`).
- Large reference assets (genome access BEDs, refFlat, gene-info TSVs) are user-supplied or downloaded, not bundled in the wheel, keeping installs small and free of build-specific or licensing baggage.

Introduce `package-data` only for a concrete need, and pair it with a CI test that installs the built wheel and loads the resource, so a missing data file fails loudly in CI rather than at a user's runtime.

## Code Style & Conventions

### Modern Python Style
- **Type annotations** use PEP 604 union syntax: `X | Y` and `X | None`, not `Union[X, Y]` or `Optional[X]`
- **Match/case** (PEP 634) is used for dispatch on string literals where it improves clarity
- **`removeprefix()`/`removesuffix()`** (PEP 616) preferred over manual slicing for prefix/suffix removal
- **Dict `|=`** (PEP 584) preferred over `.update()` for merging dict literals
- All `zip()` calls use explicit `strict=True` or `strict=False` (PEP 618)

### Imports in Tests
- `test/test_commands.py` and `test/test_cnvlib.py` each have a top-level `from cnvlib import (...)` block that serves as both a smoke test and the shared import set. Add new `cnvlib` submodule imports there rather than as local imports inside individual test methods.

### Variable Naming
- The codebase uses `bam_fname` or `sample_fname` for file paths that can be either BAM or bedGraph files
- Parameter names in function signatures often use generic terms (e.g., `bam_fname`) even when they accept multiple formats

## Design

The analytical methods implemented in CNVkit are described in the publication:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873

When implementing or modifying analytical methods, look up the primary literature to understand the underlying algorithms. Use [Google Scholar](https://scholar.google.com/) and [Europe PMC](https://europepmc.org/) to find and read the original papers for methods referenced in the code (e.g. segmentation algorithms, statistical tests, normalization approaches).
