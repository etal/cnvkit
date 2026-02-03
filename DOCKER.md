# Docker Setup for CNVkit

CNVkit provides Docker images for portable execution in pipelines and reproducible analysis.

## Using Pre-built Images

### Pull from Docker Hub

```bash
# Latest stable release
docker pull etal/cnvkit:latest

# Development version (built from master branch)
docker pull etal/cnvkit:devel

# Specific version
docker pull etal/cnvkit:0.9.13
```

### Run CNVkit in Docker

```bash
# Run CNVkit commands with your data mounted
docker run -v /path/to/data:/data etal/cnvkit:latest \
    cnvkit.py batch /data/samples/*.bam \
    -r /data/reference.cnn \
    -d /data/output

# Interactive shell
docker run -it -v /path/to/data:/data etal/cnvkit:latest bash

# Show version
docker run --rm etal/cnvkit:latest cnvkit.py version
```

## For Pipeline Users (Nextflow, WDL)

### Nextflow

```nextflow
process cnvkit_batch {
    container 'etal/cnvkit:latest'

    input:
    path bam_files
    path reference

    output:
    path "output/*"

    script:
    """
    cnvkit.py batch ${bam_files} -r ${reference} -d output
    """
}
```

### WDL

```wdl
task cnvkit_batch {
    input {
        Array[File] bam_files
        File reference
    }

    command {
        cnvkit.py batch ${sep=' ' bam_files} -r ${reference} -d output
    }

    runtime {
        docker: "etal/cnvkit:latest"
    }

    output {
        Array[File] results = glob("output/*")
    }
}
```

## Image Versioning

CNVkit Docker images follow this tagging strategy:

- **`latest`** - Most recent stable release (built from git tags)
- **`devel`** - Development version (built from master branch on each commit)
- **`X.Y.Z`** - Specific version (e.g., `0.9.13`)
- **`X.Y`** - Latest patch for minor version (e.g., `0.9` → `0.9.13`)

### Recommendation for Production

Use specific version tags for reproducibility:
```bash
docker pull etal/cnvkit:0.9.13
```

Avoid `latest` in production pipelines to ensure consistent results.

## Building Images Locally

### Build from Current Source

```bash
# Build development image
make docker

# Test it
make docker-test
```

### Build Release Image

```bash
# Build specific version from PyPI
make release_version=0.9.13 docker-release

# Push to Docker Hub (requires authentication)
make release_version=0.9.13 docker-push
```

## Automated Builds (GitHub Actions)

Docker images are automatically built and pushed by GitHub Actions:

- **On every commit to `master`**: Builds and pushes `etal/cnvkit:devel`
- **On every tag `vX.Y.Z`**: Builds and pushes:
  - `etal/cnvkit:X.Y.Z` (exact version)
  - `etal/cnvkit:X.Y` (minor version)
  - `etal/cnvkit:latest` (latest stable)

To trigger a release build:
```bash
git tag v0.9.14
git push origin v0.9.14
```

## Setting up Docker Hub Credentials

For maintainers who need to manually push images or set up CI:

1. Create Docker Hub access token at https://hub.docker.com/settings/security
2. Add secrets to GitHub repository settings:
   - `DOCKER_USERNAME`: Your Docker Hub username
   - `DOCKER_PASSWORD`: Your Docker Hub access token

## Image Contents

The production Docker image includes:

- **Python 3.10+** with CNVkit and all dependencies
- **R** with DNAcopy package for segmentation
- **Conda environment** with all required packages
- **Helper scripts** (e.g., `snpfilter.sh`)
- **Non-root user** (`cnvkit`) for security
- **Working directory** at `/data` for mounted volumes

## Troubleshooting

### Permission Issues

The container runs as user `cnvkit` (not root). If you encounter permission issues:

```bash
# Run as root (not recommended for production)
docker run --user root -v /path/to/data:/data etal/cnvkit:latest ...

# Or fix permissions on host
chmod -R a+rw /path/to/data
```

### Out of Memory

Increase Docker's memory limit:

```bash
# Run with more memory
docker run --memory=8g -v /path/to/data:/data etal/cnvkit:latest ...
```

### Version Mismatch

Verify the installed version:

```bash
docker run --rm etal/cnvkit:latest cnvkit.py version
```

## Alternative: GitHub Container Registry

While Docker Hub is the primary registry, images can also be pushed to GitHub Container Registry (ghcr.io) as a backup.

To enable, uncomment the GHCR sections in `.github/workflows/docker-build.yml` and update `DOCKER_IMAGE` to:
```yaml
DOCKER_IMAGE: ghcr.io/${{ github.repository_owner }}/cnvkit
```

## Development Container

For local development in VS Code, use the DevContainer instead:

```bash
# Build DevContainer image
make docker-dev
```

Or simply open the project in VS Code and select "Reopen in Container".

See `.devcontainer/` for configuration details.
