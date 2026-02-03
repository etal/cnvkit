# Contributing to CNVkit

Thank you for your interest in contributing to CNVkit! This document provides guidelines and information for contributors.

## Quick Links

- **Report bugs**: [GitHub Issues](https://github.com/etal/cnvkit/issues)
- **Ask questions**: [Biostars CNVkit forum](https://www.biostars.org/t/CNVkit/)
- **Documentation**: [ReadTheDocs](https://cnvkit.readthedocs.io)

## Development Setup

### Prerequisites

- Python 3.10 or later
- R with DNAcopy package (for segmentation)
- Git

### Getting Started

1. **Fork and clone the repository**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/cnvkit.git
   cd cnvkit
   ```

2. **Set up development environment** (choose one):

   **Option A: Conda (recommended)**
   ```bash
   conda env create -f environment-dev.yml
   conda activate cnvkit-dev
   pip install -e '.[test]'
   ```

   **Option B: pip**
   ```bash
   pip install -e '.[test]'
   # You'll need to install R and DNAcopy separately
   ```

   **Option C: VS Code DevContainer**
   - Open the project in VS Code
   - Select "Reopen in Container" when prompted
   - Everything is pre-configured!

3. **Install pre-commit hooks** (recommended):
   ```bash
   pre-commit install
   ```

   This automatically runs code quality checks before each commit. See `PRE-COMMIT-SETUP.md` for details.

### Making Changes

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the code style guidelines below

3. **Run tests**:
   ```bash
   # Quick tests during development
   pytest test/

   # Or use Makefile shortcuts
   make test          # Run pytest
   make lint          # Check code style
   make format        # Auto-format code
   ```

4. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Brief description of changes"
   # Pre-commit hooks will run automatically
   ```

5. **Push and create a pull request**:
   ```bash
   git push origin feature/your-feature-name
   ```
   Then open a PR on GitHub.

## Code Style and Conventions

### Python Code Style

CNVkit uses **Ruff** for linting and formatting:

- **Format**: Black-compatible (88 character line length)
- **Style**: PEP 8 with project-specific exceptions (see `pyproject.toml`)
- **Type hints**: Encouraged but not required
- **Docstrings**: Use for public APIs

Run formatting/linting:
```bash
make format    # Auto-format code
make lint      # Check for issues
```

Pre-commit hooks enforce these automatically.

### Key Conventions

- **Python version**: Code must work with Python 3.10+
- **Testing**: Use pytest, not unittest
- **Subprocess calls**: Allowed for bioinformatics tools (samtools, R, etc.)
- **Variable naming**: Use descriptive names (e.g., `bam_fname` for file paths)
- **Command architecture**: Each command has:
  - `_cmd_*()` function for CLI parsing
  - `do_*()` function for core logic/API

See `CLAUDE.md` for detailed architecture information.

## Testing

### Running Tests

```bash
# Fast: Run tests directly with pytest
pytest test/                    # All tests
pytest test/test_cnvlib.py     # Specific file
pytest -v -k "batch"           # Tests matching pattern

# Comprehensive: Run full test suite
make test-all                  # Run tox (all Python versions)
cd test && make               # Integration tests with real data
cd test && make mini          # Quick integration tests
```

### Writing Tests

- Add unit tests for new functionality in `test/`
- Use descriptive test names: `test_batch_with_bedgraph_input()`
- Test edge cases and error conditions
- Integration tests go in `test/Makefile`

### Test Coverage

We aim for high test coverage. Check coverage locally:
```bash
pytest --cov=cnvlib --cov=skgenome test/
# Or via tox
tox -e coverage
```

## Documentation

### Code Documentation

- Add docstrings to public functions and classes
- Use RST format for consistency with Sphinx
- Explain the "why", not just the "what"

### User Documentation

User-facing documentation lives in the `doc/` directory and is built with Sphinx:

```bash
make docs         # Build documentation
make docs-serve   # Build and serve locally at http://localhost:8000
```

Documentation is automatically published to [ReadTheDocs](https://cnvkit.readthedocs.io).

## Pull Request Process

1. **Ensure tests pass**: All tests must pass before merging
2. **Update documentation**: If you've added/changed features
3. **Add changelog entry**: For user-facing changes
4. **Describe your changes**: Clear PR description with context
5. **Respond to review**: Address feedback promptly

### PR Checklist

- [ ] Tests pass locally (`pytest test/`)
- [ ] Code is formatted (`make format`)
- [ ] Pre-commit hooks pass
- [ ] Documentation updated (if needed)
- [ ] Changelog entry added (if user-facing)

## Common Development Tasks

### Building a Distribution

```bash
make build    # Creates wheel and sdist in dist/
```

### Docker

```bash
make docker        # Build development Docker image
make docker-test   # Test the Docker image
```

See `DOCKER.md` for production Docker workflows.

### Security Scanning

```bash
make security     # Run safety + bandit
tox -e security   # Via tox
```

## Getting Help

- **Questions about contributing**: Open a [GitHub Discussion](https://github.com/etal/cnvkit/discussions)
- **Questions about CNVkit usage**: Ask on [Biostars](https://www.biostars.org/t/CNVkit/)
- **Bug reports**: [GitHub Issues](https://github.com/etal/cnvkit/issues)
- **Development questions**: Mention @claude in an issue for AI assistance

## Code of Conduct

- Be respectful and constructive
- Focus on what is best for the community
- Show empathy towards other community members

This is scientific software used in medical research. Please ensure:
- Correctness over cleverness
- Clear code over compact code
- Thorough testing for critical functionality

## License

By contributing to CNVkit, you agree that your contributions will be licensed under the Apache License 2.0.

## Recognition

Contributors are recognized in:
- Git commit history
- GitHub contributors page
- Release notes (for significant contributions)

Thank you for making CNVkit better! 🧬
