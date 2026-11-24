# Contributing to CNVkit

Thank you for your interest in contributing to CNVkit!

## Quick links

- **Report bugs / request features**: [GitHub Issues](https://github.com/etal/cnvkit/issues)
- **Ask usage questions**: [Biostars CNVkit forum](https://www.biostars.org/t/CNVkit/)
- **User documentation**: [ReadTheDocs](https://cnvkit.readthedocs.io)

## Getting started

The development environment is a [pixi](https://pixi.sh) workspace, declared in `pyproject.toml` and pinned by the committed lockfile `pixi.lock`.
Install pixi, then:

```bash
git clone https://github.com/YOUR_USERNAME/cnvkit.git
cd cnvkit
pixi install                     # Python, R with DNAcopy, and CNVkit in editable mode
pixi run pre-commit install      # run quality checks on commit
pixi run pytest test/            # confirm the suite passes
```

`pixi shell` starts a shell inside the environment, after which the `pixi run` prefix is unnecessary.

## Full developer guide

Environment setup, code style and conventions, testing, the pull request
process, and the maintainer workflows for building documentation and Docker
images are all documented in the **[Development and contributing
guide](https://cnvkit.readthedocs.io/en/latest/development.html)** (source:
[`doc/development.rst`](doc/development.rst)).

## License

By contributing to CNVkit, you agree that your contributions will be licensed
under the Apache License 2.0.
