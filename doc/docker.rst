Running CNVkit with Docker
==========================

CNVkit provides Docker images for portable execution in pipelines and
reproducible analysis.


Using pre-built images
----------------------

Pull the image from `Docker Hub <https://hub.docker.com/r/etal/cnvkit/>`_::

    # Latest stable release
    docker pull etal/cnvkit:latest

    # Development version (built from the master branch)
    docker pull etal/cnvkit:devel

    # A specific version
    docker pull etal/cnvkit:0.9.13

Run CNVkit commands with your data mounted into the container::

    docker run -v /path/to/data:/data etal/cnvkit:latest \
        cnvkit.py batch /data/samples/*.bam \
        -r /data/reference.cnn \
        -d /data/output

Open an interactive shell, or just print the version::

    docker run -it -v /path/to/data:/data etal/cnvkit:latest bash
    docker run --rm etal/cnvkit:latest cnvkit.py version


Use in workflow pipelines
-------------------------

The image works directly as the container for Nextflow, WDL, and similar
workflow managers.

Nextflow::

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

WDL::

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

If you would prefer a fully managed workflow, several open-source CNVkit
pipelines are available on `Dockstore
<https://dockstore.org/search?entryType=workflows&search=cnvkit>`_.


Use on Galaxy
-------------

The Intergalactic Utilities Commission (IUC) maintains a complete suite of
CNVkit tools for the `Galaxy <https://galaxyproject.org/>`_
platform, with a separate wrapper for each CNVkit subcommand (``batch``,
``access``, ``antitarget``, ``autobin``, ``coverage``, ``reference``, ``fix``,
``segment``, ``call``, ``scatter``, ``diagram``, ``heatmap``, and others).
The tools run on the public `usegalaxy.eu
<https://usegalaxy.eu/?tool_id=cnvkit_batch>`_ server and can be installed into
any Galaxy instance from the Tool Shed. Their dependencies, including CNVkit
itself, are resolved automatically through Bioconda, so no manual installation
is required. Source and issue tracking for these wrappers live in the
`tools-iuc repository
<https://github.com/galaxyproject/tools-iuc/tree/main/tools/cnvkit>`_.


Image versioning
----------------

CNVkit Docker images use the following tags:

- ``latest`` -- the most recent stable release (built from git tags)
- ``devel`` -- the development version (built from the ``master`` branch on
  each commit)
- ``X.Y.Z`` -- a specific release (e.g. ``0.9.13``)
- ``X.Y`` -- the latest patch for a minor version (e.g. ``0.9`` resolves to
  ``0.9.13``)

For reproducible analyses, pin a specific version tag rather than ``latest``::

    docker pull etal/cnvkit:0.9.13


Image contents
--------------

The production image includes:

- Python 3.11+ with CNVkit and all dependencies
- R with the DNAcopy package for segmentation
- A Conda environment with all required packages
- Helper scripts (e.g. ``snpfilter.sh``)
- A non-root user (``cnvkit``) for security
- A working directory at ``/data`` for mounted volumes


Troubleshooting
---------------

**Permission issues.** The container runs as the unprivileged user ``cnvkit``
rather than root. If you hit permission errors on mounted data, either run as
root (not recommended for production) or relax permissions on the host::

    docker run --user root -v /path/to/data:/data etal/cnvkit:latest ...
    chmod -R a+rw /path/to/data

**Out of memory.** Raise Docker's memory limit::

    docker run --memory=8g -v /path/to/data:/data etal/cnvkit:latest ...

**Version mismatch.** Confirm the version inside the image::

    docker run --rm etal/cnvkit:latest cnvkit.py version


Building images locally and pushing releases is covered in the
:doc:`development guide <development>`.
