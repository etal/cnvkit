# CNVkit Maintenance Commands
# This Makefile provides shortcuts for common development and release tasks

PYTHON := python3
PIP := pip3
TWINE := twine
DOCKER := docker

# Version for releases (use: make release_version=0.9.13 upload-pypi)
release_version ?= latest

# Docker image name
DOCKER_IMAGE := etal/cnvkit

.PHONY: help
help:
	@echo "CNVkit Maintenance Commands"
	@echo ""
	@echo "Development:"
	@echo "  install-dev    Install CNVkit in editable mode with dev dependencies"
	@echo "  install-test   Install CNVkit in editable mode with test dependencies"
	@echo "  test           Run pytest tests"
	@echo "  test-all       Run full tox test suite"
	@echo "  lint           Run ruff linting"
	@echo "  format         Run ruff formatting"
	@echo "  security       Run security checks (safety + bandit)"
	@echo "  pre-commit     Install and run pre-commit hooks"
	@echo ""
	@echo "Building & Distribution:"
	@echo "  build          Build wheel and sdist using modern build tools"
	@echo "  clean          Remove build artifacts"
	@echo ""
	@echo "Docker:"
	@echo "  docker           Build 'devel' image from current source"
	@echo "  docker-release   Build versioned release image (set release_version=X.Y.Z)"
	@echo "  docker-push      Push versioned images to Docker Hub"
	@echo "  docker-dev       Build DevContainer image for VS Code"
	@echo "  docker-test      Test the built Docker image"
	@echo "  Note: GHA automatically builds 'devel' on master, 'latest' on tags"
	@echo ""
	@echo "Release:"
	@echo "  upload-pypi    Upload to PyPI (set release_version=X.Y.Z)"
	@echo ""

# =============================================================================
# Development
# =============================================================================

.PHONY: install-dev
install-dev:
	$(PIP) install -e '.[test]'
	$(PIP) install -r requirements/dev.txt

.PHONY: install-test
install-test:
	$(PIP) install -e '.[test]'

.PHONY: test
test:
	pytest -v test/

.PHONY: test-all
test-all:
	tox

.PHONY: lint
lint:
	ruff check cnvlib skgenome

.PHONY: format
format:
	ruff format cnvlib skgenome
	ruff check --fix cnvlib skgenome

.PHONY: security
security:
	safety check
	bandit -r cnvlib skgenome -c pyproject.toml

.PHONY: pre-commit
pre-commit:
	pre-commit install
	pre-commit run --all-files

# =============================================================================
# Building & Distribution
# =============================================================================

.PHONY: build
build: clean
	$(PYTHON) -m build

.PHONY: clean
clean:
	rm -rf build dist *.egg-info
	rm -rf .tox .pytest_cache .coverage htmlcov coverage.xml
	rm -rf cnvlib/__pycache__ skgenome/__pycache__
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete

# =============================================================================
# Docker
# =============================================================================

.PHONY: docker
docker:
	@echo "Building production Docker image (from source)..."
	$(DOCKER) build -t $(DOCKER_IMAGE):devel -f Dockerfile .

.PHONY: docker-release
docker-release:
	@if [ "$(release_version)" = "latest" ]; then \
		echo "ERROR: Please specify release_version (e.g., make release_version=0.9.13 docker-release)"; \
		exit 1; \
	fi
	@echo "Building production Docker image for version $(release_version)..."
	$(DOCKER) build --build-arg CNVKIT_VERSION=$(release_version) \
		-t $(DOCKER_IMAGE):$(release_version) \
		-t $(DOCKER_IMAGE):latest \
		-f Dockerfile .

.PHONY: docker-push
docker-push:
	@if [ "$(release_version)" = "latest" ]; then \
		echo "ERROR: Please specify release_version (e.g., make release_version=0.9.13 docker-push)"; \
		exit 1; \
	fi
	@echo "Pushing Docker images for version $(release_version)..."
	$(DOCKER) push $(DOCKER_IMAGE):$(release_version)
	$(DOCKER) push $(DOCKER_IMAGE):latest

.PHONY: docker-dev
docker-dev:
	@echo "Building development DevContainer image..."
	$(DOCKER) build -t $(DOCKER_IMAGE):devcontainer -f .devcontainer/Dockerfile .

.PHONY: docker-test
docker-test:
	@echo "Testing Docker image..."
	$(DOCKER) run --rm $(DOCKER_IMAGE):devel cnvkit.py version
	$(DOCKER) run --rm $(DOCKER_IMAGE):devel cnvkit.py --help

# =============================================================================
# Release
# =============================================================================

.PHONY: upload-pypi
upload-pypi: build
	@if [ "$(release_version)" = "latest" ]; then \
		echo "ERROR: Please specify release_version (e.g., make release_version=0.9.13 upload-pypi)"; \
		exit 1; \
	fi
	@echo "Uploading version $(release_version) to PyPI..."
	$(TWINE) upload dist/*

# =============================================================================
# Documentation
# =============================================================================

.PHONY: docs
docs:
	tox -e doc

.PHONY: docs-serve
docs-serve: docs
	@echo "Serving documentation at http://localhost:8000"
	cd .tox/doc/html && $(PYTHON) -m http.server
