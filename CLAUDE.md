# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CNVkit is a command-line toolkit and Python library for detecting copy number variants and alterations genome-wide from high-throughput sequencing data. It provides both a CLI interface and Python API for genomic analysis workflows.

**Supported Python versions:** 3.10+ (tested on 3.10-3.14)
**Minimum versions aligned with Ubuntu 25.04 (Plucky)**

## Development Workflow

- **Bug fixes and new features**: Write a failing test first, then implement.
- **Edge cases**: Before finishing, verify behavior for empty inputs, NaN/missing values, and single-element arrays.
- **NaN weight safety**: `np.average(values, weights=w)` and `numpy.sum(w)` propagate NaN; use `np.nansum` for weight sums and filter `~np.isnan(wt)` before `np.average`. Note that `pandas.Series.sum()` skips NaN by default — prefer explicit `np.nansum` for clarity.
- **User-facing changes**: Update the relevant docs in `doc/*.rst`.
- **Clinical impact**: When reviewing changes, consider whether the changeset alters numerical output or output file formats (`.cnr`, `.cns`, `.cnn`, SEG, VCF). Flag any such changes explicitly, as downstream clinical pipelines may depend on exact output stability.

## Development Commands

Use 'bd' for task tracking.

### Development Environment

**Conda (recommended):**
```bash
conda env create -f environment-dev.yml   # All deps, R, testing tools
conda activate cnvkit
```

**Note:** The conda env must be activated before running pytest, mypy, or other dev tools -- they are not installed globally. The conda env includes R with DNAcopy for segmentation. Use `conda activate cnvkit && <command>` in scripts; `conda run` is unreliable here.

### Testing

Run tests directly with pytest (not tox):
```bash
pytest test/                           # Run all tests
pytest test/test_cnvlib.py            # Run specific test file
pytest test/test_commands.py::CommandTests::test_batch  # Run specific test
pytest test/test_commands.py -k "batch or coverage"  # Run tests matching patterns
```

**Comprehensive testing:**
- `tox` - Full matrix: Python 3.10-3.14, linting, security, coverage, docs
- `cd test/ && make mini` - Integration tests with real genomic data (used in CI)

### Type Checking

```bash
mypy                              # Check both cnvlib and skgenome
mypy cnvlib/batch.py              # Check a single file
```

**mypy configuration** (`pyproject.toml [tool.mypy]`):
- `check_untyped_defs = true`, `warn_unreachable = true`, `warn_return_any = true`
- `enable_error_code = ["ignore-without-code"]` - All `# type: ignore` must include error codes

**Common type patterns in this codebase:**
- `tabio.read()` returns union types -- use `# type: ignore[return-value]` at call sites
- Pandas operations often return `Any` -- use `# type: ignore[no-any-return]`
- Closures don't narrow types in mypy -- use `# type: ignore[index]` or local asserts
- Parameters typed `param: None = None` cause unreachable blocks -- use `Optional[Type] = None` instead
- Generator functions must use `-> Generator[YieldType, SendType, ReturnType]`, not `-> Iterator` or `-> None`
- When a variable changes type (e.g. `str` → `list[int]`), rename it to avoid shadowing (e.g. `copies` → `copy_strs`)
- Use `assert x is not None` to narrow `Optional` types when control flow guarantees non-None
- `numpy.bool_` is not assignable to `bool` in mypy -- use `# type: ignore[assignment]` or widen parameters to `bool | bool_ | None`

### Code Quality & Security

Pre-commit hooks run automatically on `git commit` (ruff, bandit, whitespace, YAML/TOML checks). Install with `pre-commit install`.

**Manual commands:**
- `make lint` / `make format` - Ruff linting and formatting
- `make security` - Safety + Bandit scans
- `tox -e lint` / `tox -e security` / `tox -e coverage` / `tox -e doc`
- `make help` - Show all Makefile targets

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

GenomicArray uses a `TypeVar` pattern (`_GA = TypeVar("_GA", bound="GenomicArray")`) so that methods like `.copy()`, `.concat()`, and `.as_dataframe()` preserve the subclass type through type checking.

### File Formats
- `.cnn` - Coverage/reference data
- `.cnr` - Copy number ratio data
- `.cns` - Segmented copy number data

## Dependencies

Core dependencies are declared in `requirements/core.txt`; `min.txt` pins exact minimums for compatibility testing.

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

## Serena (MCP LSP integration)

A Serena MCP server provides LSP-backed code intelligence tools (configured via `claude mcp add` with `--context claude-code`, which exposes only LSP tools, not file I/O or shell).

**When to use Serena vs. built-in tools:**
- **Exploring unfamiliar code** — use `get_symbols_overview` to see a file's structure without reading the whole file, then `find_symbol` with `include_body=True` to read only the methods you need
- **Tracing call graphs** — use `find_referencing_symbols` to find all callers/users of a symbol across the codebase
- **Refactoring** — use `replace_symbol_body` / `insert_after_symbol` for precise symbolic edits
- **Simple lookups** — use Grep/Glob for known string patterns; Serena is better for semantic queries (e.g. "all methods of CopyNumArray" or "all callers of `do_segmentation`")

**Key tools:**
- `find_symbol` - Find by name path (e.g. `CopyNumArray/squash_genes`); use `depth=1` to list methods, `include_body=True` to read implementations
- `get_symbols_overview` - List all top-level symbols in a file
- `find_referencing_symbols` - Find all references to a symbol with surrounding code context
- `search_for_pattern` - Regex search with file filtering (for non-symbol searches)
- `replace_symbol_body` / `insert_after_symbol` / `insert_before_symbol` - Symbolic edits

## Design

The analytical methods implemented in CNVkit are described in the publication:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873

When implementing or modifying analytical methods, look up the primary literature to understand the underlying algorithms. Use [Google Scholar](https://scholar.google.com/) and [Europe PMC](https://europepmc.org/) to find and read the original papers for methods referenced in the code (e.g. segmentation algorithms, statistical tests, normalization approaches).
