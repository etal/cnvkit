Development and contributing
============================

Thank you for your interest in contributing to CNVkit! This page is the
detailed developer guide: environment setup, code style, testing, the pull
request process, and the maintainer workflows for building documentation and
Docker images.

CNVkit is scientific software used in medical research, so we value:

- Correctness over cleverness
- Clear code over compact code
- Thorough testing for critical functionality

For day-to-day architecture and coding conventions, see ``CLAUDE.md`` in the
repository root.


Development setup
-----------------

Prerequisites:

- Python 3.11 or later
- R with the DNAcopy package (for the default CBS segmentation)
- Git

Set up a development environment using one of the following options.

**Option A: Conda (recommended).** This brings in R and DNAcopy as well::

    conda env create -f environment-dev.yml
    conda activate cnvkit
    pip install -e '.[test]'

**Option B: pip.** You will need to install R and DNAcopy separately::

    pip install -e '.[test]'

**Option C: VS Code DevContainer.** Open the project in VS Code and choose
"Reopen in Container"; everything is pre-configured.

Finally, install the pre-commit hooks so quality checks run on every commit::

    pre-commit install


Pre-commit hooks
----------------

The project uses `pre-commit <https://pre-commit.com/>`_ to run code-quality
checks automatically before each commit. The configuration lives in
``.pre-commit-config.yaml``, and the tools it runs are configured in
``pyproject.toml`` (the ``[tool.ruff]`` and ``[tool.bandit]`` sections).

What gets checked
~~~~~~~~~~~~~~~~~

On every ``git commit``:

1. **Ruff linting** -- style (PEP 8), common bugs, unused imports/variables,
   and modern Python idioms (pyupgrade)
2. **Ruff formatting** -- Black-compatible auto-formatting
3. **Bandit** -- scans for common security issues
4. **File hygiene** -- trailing whitespace, end-of-file newline, YAML/TOML
   syntax, merge-conflict markers, and large files (>1 MB)

When a hook fails, the commit is blocked. Formatting and whitespace issues are
auto-fixed in place -- review the changes with ``git diff``, ``git add`` them,
and commit again.

Running hooks manually
~~~~~~~~~~~~~~~~~~~~~~~

::

    pre-commit run                 # On staged files
    pre-commit run --all-files     # On the whole tree
    pre-commit run ruff --all-files  # A single hook
    pre-commit autoupdate          # Update hooks to their latest versions

The root ``Makefile`` also wraps the underlying tools::

    make pre-commit    # Install hooks and run on all files
    make lint          # Run ruff
    make format        # Auto-format code
    make security      # Run pip-audit + bandit

Skipping hooks (sparingly)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    git commit --no-verify    # Skip all hooks (not recommended)
    SKIP=bandit git commit    # Skip a specific hook

If a hook environment misbehaves, clear the pre-commit cache::

    pre-commit clean
    pre-commit gc


Code style and conventions
--------------------------

CNVkit uses **Ruff** for both linting and formatting:

- **Format:** Black-compatible, 88-character line length
- **Style:** PEP 8 with project-specific exceptions (see ``pyproject.toml``)
- **Type hints:** encouraged; checked with mypy (``mypy`` checks ``cnvlib`` and
  ``skgenome``)
- **Docstrings:** use for public APIs, in RST format for Sphinx

Key conventions:

- Code must work with Python 3.11+
- Use pytest, not unittest
- Subprocess calls are allowed for bioinformatics tools (samtools, R, etc.)
- Use descriptive variable names (e.g. ``bam_fname`` for file paths)
- Each command has a ``_cmd_*()`` function for CLI parsing and a ``do_*()``
  function for the core logic / Python API

See ``CLAUDE.md`` for the full architecture notes and modern-Python style rules.


Testing
-------

Run the test suite directly with pytest::

    pytest test/                   # All tests
    pytest test/test_cnvlib.py     # A specific file
    pytest -v -k "batch"           # Tests matching a pattern

Or use the Makefile and tox for broader coverage::

    make test          # pytest -n auto
    make test-all      # tox: full matrix across Python versions
    cd test && make mini   # Integration tests with real genomic data

When writing tests:

- Add unit tests for new functionality under ``test/``
- Write a failing test first, then implement the fix or feature
- Use descriptive names, e.g. ``test_batch_with_bedgraph_input()``
- Cover edge cases: empty inputs, NaN/missing values, single-element arrays
- Integration tests go in ``test/Makefile``

Check coverage locally::

    pytest --cov=cnvlib --cov=skgenome test/
    tox -e coverage


Documentation
-------------

User-facing documentation lives in this ``doc/`` directory and is built with
Sphinx, then published to `ReadTheDocs <https://cnvkit.readthedocs.io>`_::

    make docs          # Build the HTML docs (via tox -e doc)
    make docs-serve    # Build and serve locally at http://localhost:8000

When you add or change a user-facing feature, update the relevant ``doc/*.rst``
page in the same pull request.


Pull request process
--------------------

1. Create a feature branch (``git checkout -b feature/your-feature-name``).
2. Make your changes following the conventions above.
3. Ensure the tests pass and the code is formatted.
4. Update documentation for any user-facing change.
5. Push your branch and open a pull request on GitHub with a clear description.

Before opening a PR, confirm:

- Tests pass locally (``pytest test/``)
- Code is formatted (``make format``)
- Pre-commit hooks pass
- Documentation is updated (if needed)
- A changelog entry is added (for user-facing changes)


Building distributions
----------------------

Build a wheel and source distribution into ``dist/``::

    make build


Releasing to PyPI
-----------------

The full, narrative release procedure lives in the `Release procedure
<https://github.com/etal/cnvkit/wiki/Release-procedure>`_ wiki page. The
repository enforces one invariant mechanically: **the ``master`` branch always
carries a development version** (``X.Y.Z.devN``), and a plain ``X.Y.Z`` version
exists only transiently on the tagged release commit. The ``Version Guard``
GitHub Actions workflow fails any master build or pull request whose
``cnvlib/_version.py`` is left on a plain release version, which prevents the
post-release ``.dev0`` bump from being silently skipped. Check it locally with::

    make check-version

To cut a release, once the notes are written:

1. Set ``cnvlib/_version.py`` to the plain release version ``X.Y.Z`` (drop the
   ``.devN`` suffix) and commit.
2. Tag the commit with the summarized changelog as its message and push::

       git tag vX.Y.Z -aF release-notes.md
       git push --tags

   GitHub Actions builds and publishes the Docker image for the tag.
3. Build and upload the Python distributions to PyPI::

       make build
       make release_version=X.Y.Z upload-pypi

4. Bump ``cnvlib/_version.py`` to the next development version
   ``X.Y.(Z+1).dev0``, commit, and push. (``Version Guard`` enforces that this
   step is not forgotten.)


Building and releasing Docker images
------------------------------------

End-user instructions for *running* the published images are in
:doc:`docker`. The workflows below are for maintainers.

Build and smoke-test images locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    make docker        # Build the production image (etal/cnvkit:devel) from source
    make docker-test   # Run cnvkit.py version inside the image
    make docker-dev    # Build the DevContainer image

Build and push a release
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    make release_version=0.9.14 docker-release   # Build a tagged image from PyPI
    make release_version=0.9.14 docker-push      # Push to Docker Hub (needs auth)

Automated builds (GitHub Actions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Docker images are built and pushed automatically:

- On every commit to ``master``: ``etal/cnvkit:devel``
- On every tag ``vX.Y.Z``: ``etal/cnvkit:X.Y.Z``, ``etal/cnvkit:X.Y``, and
  ``etal/cnvkit:latest``

To trigger a release build, push a version tag::

    git tag v0.9.14
    git push origin v0.9.14

Docker Hub credentials
~~~~~~~~~~~~~~~~~~~~~~~

For maintainers who push images or set up CI:

1. Create a Docker Hub access token at
   https://hub.docker.com/settings/security
2. Add these secrets to the GitHub repository settings:

   - ``DOCKER_USERNAME`` -- your Docker Hub username
   - ``DOCKER_PASSWORD`` -- your Docker Hub access token

GitHub Container Registry (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Docker Hub is the primary registry, but images can also be pushed to
``ghcr.io`` as a backup. To enable it, uncomment the GHCR sections in
``.github/workflows/docker-build.yml`` and set ``DOCKER_IMAGE`` to
``ghcr.io/${{ github.repository_owner }}/cnvkit``.


Getting help
------------

- Questions about contributing: open a
  `GitHub Discussion <https://github.com/etal/cnvkit/discussions>`_
- Questions about CNVkit usage: ask on
  `Biostars <https://www.biostars.org/t/CNVkit/>`_
- Bug reports: `GitHub Issues <https://github.com/etal/cnvkit/issues>`_


Code of conduct
---------------

Be respectful and constructive, focus on what is best for the community, and
show empathy toward other contributors.


License
-------

By contributing to CNVkit, you agree that your contributions will be licensed
under the Apache License 2.0.
