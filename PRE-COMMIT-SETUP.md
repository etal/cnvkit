# Pre-commit Setup Guide

This project now uses pre-commit hooks to automatically check code quality before commits.

## Quick Start

```bash
# 1. Install pre-commit (if not already installed)
brew install pre-commit
# or: pip install pre-commit

# 2. Install the git hooks
pre-commit install

# 3. (Optional) Run on all files to test
pre-commit run --all-files
```

## What Gets Checked

The pre-commit hooks will automatically run on every `git commit`:

1. **Ruff Linting** - Fast Python linter checking for:
   - Code style issues (PEP 8)
   - Common bugs and anti-patterns
   - Unused imports and variables
   - Modern Python idioms (pyupgrade)

2. **Ruff Formatting** - Automatic code formatting (Black-compatible)

3. **Security Checks** - Bandit scans for common security issues

4. **File Quality**:
   - Remove trailing whitespace
   - Ensure files end with newline
   - Check YAML/TOML syntax
   - Detect merge conflicts
   - Check for large files (>1MB)

## Usage

### Automatic (Recommended)
Once installed, hooks run automatically on `git commit`. If any hook fails:
- The commit is blocked
- Files are auto-fixed when possible (formatting, whitespace)
- Review the changes with `git diff`
- Stage fixes with `git add`
- Commit again

### Manual
```bash
# Run on staged files
pre-commit run

# Run on all files
pre-commit run --all-files

# Run specific hook
pre-commit run ruff --all-files

# Update hooks to latest versions
pre-commit autoupdate
```

### Skip Hooks (Use Sparingly)
```bash
# Skip all hooks for a commit (not recommended)
git commit --no-verify

# Skip specific hook
SKIP=bandit git commit
```

## Makefile Integration

The root `Makefile` provides shortcuts:

```bash
make pre-commit    # Install hooks and run on all files
make lint          # Run ruff manually
make format        # Auto-format code
make security      # Run safety + bandit
```

## Configuration

- **Pre-commit config**: `.pre-commit-config.yaml`
- **Ruff config**: `pyproject.toml` (see `[tool.ruff]` section)
- **Bandit config**: `pyproject.toml` (see `[tool.bandit]` section)

## Troubleshooting

### Hooks fail on test files
Test files are excluded from bandit security checks. If you're getting false positives, add exclusions to `pyproject.toml`.

### Ruff auto-fix changes code unexpectedly
Review the changes with `git diff`. Ruff follows Black formatting style. If needed, adjust settings in `pyproject.toml`.

### Large files or slow performance
Pre-commit caches hook environments in `~/.cache/pre-commit/`. Clear with:
```bash
pre-commit clean
pre-commit gc
```

### Update to latest hook versions
```bash
pre-commit autoupdate
```

## CI Integration

Pre-commit hooks complement (but don't replace) CI checks:
- **Local**: Fast feedback on common issues before push
- **CI**: Comprehensive testing across Python versions, security scans, coverage

Both use the same tools (ruff, bandit) for consistency.
