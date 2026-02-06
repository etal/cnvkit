# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CNVkit is a command-line toolkit and Python library for detecting copy number variants and alterations genome-wide from high-throughput sequencing data. It provides both a CLI interface and Python API for genomic analysis workflows.

**Supported Python versions:** 3.10+ (tested on 3.10-3.14)
**Minimum versions aligned with Ubuntu 25.04 (Plucky)**

## Development Commands

### Development Environment

**Conda (recommended):**
For development in a terminal-based environment (e.g. vim/neovim), use the conda environment file:
- `environment-dev.yml` - Complete development environment (Python 3.11, all deps, R, testing tools)

```bash
conda env create -f environment-dev.yml
conda activate cnvkit-dev
```

**Note:** The conda env must be activated before running pytest, mypy, or other dev tools -- they are not installed globally.

**Pip/Manual Setup:**
Alternatively, use pip to install CNVkit in editable mode:
```bash
pip install -e '.[test]'
```

You'll still need to install R and DNAcopy separately for segmentation (see R Dependencies below).

**VS Code with DevContainer:**
The project includes a devcontainer configuration that provides a pre-configured development environment:
- All Python dependencies pre-installed via conda
- R and DNAcopy package pre-configured
- CNVkit installed in editable mode
- Simply open the project in VS Code and select "Reopen in Container" when prompted

**Docker (for production/pipeline use):**
Pre-built Docker images are available for portable execution. See `DOCKER.md` for details:
- `docker pull etal/cnvkit:latest` - Latest stable release
- `docker pull etal/cnvkit:devel` - Development version from master
- Automatically built via GitHub Actions on every commit and release


### Testing

**For local development iteration:**
Run tests directly with pytest (not tox):
```bash
# From project root - prefix paths with test/
pytest test/                           # Run all tests
pytest test/test_cnvlib.py            # Run specific test file
pytest test/test_commands.py::CommandTests::test_batch  # Run specific test
pytest -v test/                        # Verbose output
pytest test/test_commands.py -k "batch or coverage"  # Run tests matching patterns

# From test directory - no test/ prefix needed
cd test
pytest test_commands.py                # Run specific test file
pytest test_commands.py::CommandTests::test_batch  # Run specific test
pytest test_commands.py -k bedgraph    # Run tests matching pattern
pytest test_commands.py -v -k "batch or coverage"  # Multiple patterns, verbose
pytest test_commands.py --collect-only # List available tests without running
```

**For comprehensive testing:**
- `tox` - Run all environments: tests across Python versions, linting, security, coverage, docs
- `cd test/ && make` - Run comprehensive integration tests using real genomic data
- `cd test/ && make mini` - Run minimal integration tests (used in CI)

The full tox matrix (multiple Python versions, min/max deps) runs in CI. For development iteration, use pytest directly.

### Type Checking

The project uses mypy for static type checking, configured in `pyproject.toml`:
```bash
mypy                              # Check both cnvlib and skgenome
mypy cnvlib/batch.py              # Check a single file
tox -e typecheck                  # Run via tox
```

**mypy configuration** (`pyproject.toml [tool.mypy]`):
- `check_untyped_defs = true` - Check function bodies even without annotations
- `warn_unreachable = true` - Flag dead code from type narrowing
- `warn_return_any = true` - Flag functions that return `Any` unexpectedly
- `enable_error_code = ["ignore-without-code"]` - Require specific error codes in `# type: ignore` comments

**Common type patterns in this codebase:**
- `tabio.read()` returns union types -- use `# type: ignore[return-value]` at call sites
- Pandas operations often return `Any` -- use `# type: ignore[no-any-return]`
- Closures don't narrow types in mypy -- use `# type: ignore[index]` or local asserts
- Parameters typed `param: None = None` cause unreachable blocks -- use `Optional[Type] = None` instead
- Generator functions must use `-> Generator[YieldType, SendType, ReturnType]`, not `-> Iterator` or `-> None`
- When a variable changes type (e.g. `str` → `list[int]`), rename it to avoid shadowing (e.g. `copies` → `copy_strs`)
- Use `assert x is not None` to narrow `Optional` types when control flow guarantees non-None

### Code Quality & Security

**Pre-commit hooks (recommended for contributors):**
```bash
# Install pre-commit hooks (one-time setup)
pre-commit install

# Run manually on all files
pre-commit run --all-files

# Hooks run automatically on git commit
```

The pre-commit configuration includes:
- Ruff linting and formatting (including TC003: stdlib imports in TYPE_CHECKING blocks)
- Trailing whitespace removal
- YAML/TOML validation
- Bandit security checks
- Python best practices checks (including: use `x: list[str] = []` not `x = []  # type: list[str]`)

**Manual commands:**
- `make lint` - Run ruff linting
- `make format` - Auto-format code with ruff
- `make security` - Run security scans (safety + bandit)
- `make pre-commit` - Install and run pre-commit hooks
- `tox -e lint` - Run ruff via tox
- `tox -e security` - Run security scans via tox
- `tox -e coverage` - Run tests with coverage reporting

### Documentation
- `tox -e doc` - Build Sphinx documentation
- Documentation source is in `doc/` directory

### Building/Distribution
- `make build` - Build wheel and sdist using modern build tools (python -m build)
- `make clean` - Remove build artifacts and caches
- `make help` - Show all available Makefile targets
- `tox -e build` - Build wheel package via tox

### Quick Start with Makefile
The root `Makefile` provides convenient shortcuts:
```bash
make help          # Show all available commands
make install-dev   # Install in editable mode with all dev dependencies
make test          # Run pytest on test directory
make test-all      # Run full tox suite (all Python versions)
make lint          # Run ruff linting
make format        # Auto-format with ruff
make security      # Run safety + bandit
make pre-commit    # Setup and run pre-commit hooks
make build         # Build distribution packages
make clean         # Clean all build artifacts
```

## Architecture

### Core Structure
- **`cnvlib/`** - Main Python package containing the CNVkit library
  - `cnvkit.py` - CLI entry point that routes to commands
  - `commands.py` - Command-line interface definitions and API functions
  - `core.py` - Core data structures and utilities
  - `segmentation/` - Segmentation algorithms (CBS, HMM, etc.)
  - `cli/` - Additional CLI utilities and scripts

- **`skgenome/`** - Genomic data handling library (part of CNVkit but decoupled)
  - `gary.py` - GenomicArray class for genomic interval data
  - `tabio/` - File I/O for various genomic formats

### Command Architecture
CNVkit follows a pattern where each command has:
- `_cmd_*` function in commands.py - handles CLI argument parsing and I/O
- `do_*` function - implements the core functionality as an API
- Commands are auto-discovered and registered via the argparse system

### Key Data Types
- `GenomicArray` (from skgenome) - Core data structure for genomic intervals, wraps a pandas DataFrame
- `CopyNumArray` (cnary.py) - Extends GenomicArray with log2 ratios and gene names
- `VariantArray` (vary.py) - Extends GenomicArray with variant allele data
- `.cnn` files - Coverage/reference data
- `.cnr` files - Copy number ratio data
- `.cns` files - Segmented copy number data

GenomicArray uses a `TypeVar` pattern (`_GA = TypeVar("_GA", bound="GenomicArray")`) so that methods like `.copy()`, `.concat()`, and `.as_dataframe()` preserve the subclass type through type checking.

### File Format Support
The codebase handles multiple genomic formats through `skgenome.tabio`:
- BED, GFF, VCF, SEG formats
- Picard CalculateHsMetrics output
- Custom CNVkit formats (.cnn, .cnr, .cns)
- bedGraph format (.bed.gz with tabix index) - supported as input for coverage and batch commands

## Dependencies & Requirements

### Python Dependencies
Core dependencies are managed in `requirements/` with versions aligned to Ubuntu 24.04 LTS:
- `core.txt` - Essential runtime dependencies (numpy≥1.26.4, pandas≥2.1.4, etc.)
- `tests.txt` - Testing dependencies (pytest, coverage, ruff, safety, bandit)
- `dev.txt` - Development dependencies
- `min.txt` - Exact minimum versions for compatibility testing

### R Dependencies
Segmentation algorithms require R with the DNAcopy package.

**In devcontainer:** R and DNAcopy are pre-installed.

**Manual installation:**
```r
BiocManager::install("DNAcopy")
```

**Verify R setup:**
If tests fail with R-related errors, verify the installation:
```bash
Rscript -e "library(DNAcopy)"
```

## Testing & CI Strategy

### Local Testing
- **Unit tests:** `pytest test/` - Fast tests for core functionality
- **Integration tests:** `test/Makefile` - Real genomic data workflows
- **Multi-version testing:** `tox` - Tests Python 3.10-3.14 + quality checks
- **Security scanning:** `tox -e security` - Dependency vulnerabilities + static analysis

### GitHub Actions CI
- **Unit tests:** Python 3.10-3.14 on Ubuntu + macOS (3.10, 3.12)
- **Minimum version testing:** Ensures compatibility with Ubuntu LTS packages
- **Code quality:** Ruff linting with comprehensive rule set
- **Type checking:** mypy strict mode (0 errors enforced)
- **Security:** Safety (dependencies) + Bandit (code analysis)
- **Coverage:** Pytest-cov with Codecov integration
- **Integration tests:** Real genomic data pipeline (`make mini`)

### Tox Environments
- `py3{10,11,12,13,14}` - Unit tests across Python versions
- `py311-min`, `py310-min` - Test minimum dependency versions
- `lint` - Ruff code quality checks
- `typecheck` - mypy static type checking
- `security` - Safety + Bandit security scans
- `coverage` - Coverage reporting with XML output
- `doc` - Sphinx documentation build
- `build` - Package building

## Code Style & Conventions

### Variable Naming
- The codebase uses `bam_fname` or `sample_fname` for file paths that can be either BAM or bedGraph files
- When updating help text or documentation, maintain consistency with existing patterns
- Parameter names in function signatures often use generic terms (e.g., `bam_fname`) even when they accept multiple formats

## Code Quality Tools

### Ruff Configuration
The project uses ruff for unified linting, formatting, and code quality:
- **Target:** Python 3.10+
- **Line length:** 88 characters
- **Rules:** pycodestyle, pyflakes, pyupgrade, bugbear, simplify, type-checking (TC)
- **Config:** `[tool.ruff]` section in `pyproject.toml`

### Mypy Configuration
Static type checking is configured in `pyproject.toml [tool.mypy]`:
- **Target:** Python 3.10
- **Packages:** `cnvlib`, `skgenome`
- **Strict checks:** `check_untyped_defs`, `warn_unreachable`, `warn_return_any`
- **All `# type: ignore` comments must include error codes** (e.g. `# type: ignore[return-value]`)

### Security Tools
- **Safety:** Scans dependencies for known vulnerabilities
- **Bandit:** Static analysis for common security issues in Python code

### Docker Support
- **Dockerfile** with flexible version handling for deployment
- **Build args:** `CNVKIT_VERSION` for specific version installs
- **Base:** `continuumio/miniconda3:latest` with conda environment
- **DevContainer:** `.devcontainer/` provides VS Code dev environment (alternative to conda setup)

## File Organization Tips

- Main command implementations are in individual modules (e.g., `batch.py`, `segment.py`)
- Utility functions are in `cmdutil.py`, `params.py`
- Plotting/visualization code is in `plots.py`, `diagram.py`, `heatmap.py`, `scatter.py`
- Data import/export functions are in `importers.py`, `export.py`
- The `vary.py` and `cnary.py` modules extend GenomicArray for CNVkit-specific data

## Design
The analytical methods implemented in CNVkit are described in the publication:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873
