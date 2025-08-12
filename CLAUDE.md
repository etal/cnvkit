# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CNVkit is a command-line toolkit and Python library for detecting copy number variants and alterations genome-wide from high-throughput sequencing data. It provides both a CLI interface and Python API for genomic analysis workflows.

**Supported Python versions:** 3.9+ (tested on 3.9-3.13)  
**Minimum versions aligned with Ubuntu 24.04 LTS**

## Development Commands

### Testing
- `pytest test/` - Run the test suite
- `pytest test/test_cnvlib.py::test_function` - Run a specific test
- `tox` - Run all environments: tests across Python versions, linting, security, coverage, docs
- `cd test/ && make` - Run comprehensive integration tests using real genomic data
- `cd test/ && make mini` - Run minimal integration tests (used in CI)

### Code Quality & Security
- `tox -e lint` - Run ruff linting and code quality checks
- `tox -e security` - Run security scans (safety + bandit)
- `ruff check cnvlib skgenome` - Run linting manually
- `ruff format cnvlib skgenome` - Auto-format code
- `tox -e coverage` - Run tests with coverage reporting

### Documentation
- `tox -e doc` - Build Sphinx documentation
- Documentation source is in `doc/` directory

### Building/Distribution
- `tox -e build` - Build wheel package
- `make` - Run maintenance commands (see root Makefile)

### Environment Setup
For development, install in editable mode:
```bash
# Basic installation (no HMM segmentation)
pip install -e .

# With HMM segmentation support
pip install -e .[hmm]
```

Use conda environment files for dependency management:
- `conda-env.yml` - Main conda environment with Python ≥3.9 (no HMM)
- `conda-env-hmm.yml` - HMM-enabled environment with pomegranate + PyTorch
- `base.conda-env.yml` - Base environment
- `testenv.conda-env.yml` - Testing environment

## Architecture

### Core Structure
- **`cnvlib/`** - Main Python package containing the CNVkit library
  - `cnvkit.py` - CLI entry point that routes to commands
  - `commands.py` - Command-line interface definitions and API functions
  - `core.py` - Core data structures and utilities
  - `segmentation/` - Segmentation algorithms (CBS, HMM, etc.)
  - `cli/` - Additional CLI utilities and scripts

- **`skgenome/`** - Genomic data handling library (part of CNVkit)
  - `gary.py` - GenomicArray class for genomic interval data
  - `tabio/` - File I/O for various genomic formats

### Command Architecture
CNVkit follows a pattern where each command has:
- `_cmd_*` function in commands.py - handles CLI argument parsing and I/O
- `do_*` function - implements the core functionality as an API
- Commands are auto-discovered and registered via the argparse system

### Key Data Types
- `GenomicArray` (from skgenome) - Core data structure for genomic intervals
- `.cnn` files - Coverage/reference data 
- `.cnr` files - Copy number ratio data
- `.cns` files - Segmented copy number data

### File Format Support
The codebase handles multiple genomic formats through `skgenome.tabio`:
- BED, GFF, VCF, SEG formats
- Picard CalculateHsMetrics output
- Custom CNVkit formats (.cnn, .cnr, .cns)

## Dependencies & Requirements

### Python Dependencies
Core dependencies are managed in `requirements/` with versions aligned to Ubuntu 24.04 LTS:
- `core.txt` - Essential runtime dependencies (numpy≥1.26.4, pandas≥2.1.4, etc.)
- `hmm.txt` - Optional HMM dependencies (pomegranate≥1.0.0, torch≥1.13.0)
- `tests.txt` - Testing dependencies (pytest, coverage, ruff, safety, bandit)
- `dev.txt` - Development dependencies  
- `min.txt` - Exact minimum versions for compatibility testing

**Note:** HMM segmentation (`-m hmm`) requires heavy dependencies (~500MB+) and is now optional.
Install with: `pip install cnvkit[hmm]` or `conda env create -f conda-env-hmm.yml`

### R Dependencies
Segmentation algorithms require R with the DNAcopy package:
```r
BiocManager::install("DNAcopy")
```

## Testing & CI Strategy

### Local Testing
- **Unit tests:** `pytest test/` - Fast tests for core functionality
- **Integration tests:** `test/Makefile` - Real genomic data workflows
- **Multi-version testing:** `tox` - Tests Python 3.9-3.13 + quality checks
- **HMM testing:** `tox -e py312-hmm` - Test with HMM dependencies
- **Security scanning:** `tox -e security` - Dependency vulnerabilities + static analysis

### GitHub Actions CI
- **Unit tests:** Python 3.9-3.13 on Ubuntu + macOS (3.9, 3.12)
- **HMM tests:** Python 3.9-HMM, 3.12-HMM with pomegranate + PyTorch
- **Minimum version testing:** Ensures compatibility with Ubuntu LTS packages
- **Code quality:** Ruff linting with comprehensive rule set
- **Security:** Safety (dependencies) + Bandit (code analysis)  
- **Coverage:** Pytest-cov with Codecov integration (includes HMM)
- **Integration tests:** Real genomic data pipeline (`make mini`) with HMM support

### Tox Environments
- `py3{9,10,11,12,13}` - Unit tests across Python versions (core functionality)
- `py3{12,9}-hmm` - Unit tests with HMM dependencies (pomegranate + PyTorch)
- `py311-min`, `py39-min` - Test minimum dependency versions
- `lint` - Ruff code quality checks  
- `security` - Safety + Bandit security scans
- `coverage` - Coverage reporting with XML output (includes HMM)
- `doc` - Sphinx documentation build
- `build` - Package building

## Code Quality Tools

### Ruff Configuration
The project uses ruff for unified linting, formatting, and code quality:
- **Target:** Python 3.9+
- **Line length:** 88 characters
- **Rules:** pycodestyle, pyflakes, pyupgrade, bugbear, simplify
- **Config:** `[tool.ruff]` section in `pyproject.toml`

### Security Tools
- **Safety:** Scans dependencies for known vulnerabilities
- **Bandit:** Static analysis for common security issues in Python code

### Docker Support
- **Modern Dockerfile** with flexible version handling
- **Build args:** `CNVKIT_VERSION` for specific version installs
- **Base:** `continuumio/miniconda3:latest` with conda environment

## HMM Segmentation (Optional)

CNVkit v0.9.12+ includes optional HMM segmentation using pomegranate v1.0.0 with PyTorch backend.

### Key Changes from Previous Versions
- **Performance:** Reduced max_iter from 100,000 to 100 for faster convergence
- **Dependencies:** Pomegranate moved to optional `[hmm]` extra (~500MB+ with PyTorch)
- **API:** Updated to pomegranate 1.0.0 DenseHMM API with proper dtype handling
- **Error handling:** Graceful fallback when HMM dependencies unavailable

### Available Methods
- `hmm` - 3-state flexible model (loss, neutral, gain)
- `hmm-tumor` - 5-state tumor model (deletion, loss, neutral, gain, amplification)
- `hmm-germline` - 3-state fixed-means germline model

### Implementation Notes
- Located in `cnvlib/segmentation/hmm.py` with optional import handling
- Uses float32 tensors for PyTorch compatibility
- Constants defined for performance tuning (DEFAULT_MAX_ITER, DEFAULT_TOLERANCE)
- Tests conditionally skip when dependencies unavailable

## File Organization Tips

- Main command implementations are in individual modules (e.g., `batch.py`, `segment.py`)
- Utility functions are in `cmdutil.py`, `params.py`
- Plotting/visualization code is in `plots.py`, `diagram.py`, `heatmap.py`, `scatter.py`
- Data import/export functions are in `importers.py`, `export.py`
- The `vary.py` and `cnary.py` modules extend GenomicArray for CNVkit-specific data