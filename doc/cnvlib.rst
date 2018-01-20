cnvlib package
==============

Module ``cnvlib`` contents
--------------------------

.. automodule:: cnvlib
    :members:
    :undoc-members:
    :show-inheritance:

The one function exposed at the top level, `read`, loads a file in CNVkit's
BED-like tabular format and returns a CopyNumArray instance. For your own
scripting, you can usually accomplish what you need using just the CopyNumArray
and GenomicArray methods available on this returned object (see
:ref:`api_classes`).

To load other file formats, see :ref:`api_tabio`. To run functions equivalent to
CNVkit commands within Python, see :ref:`api_commands`.


.. _api_classes:

Core classes
------------

The core objects used throughout CNVkit. The base class is GenomicArray from
:doc:`skgenome`. All of these classes wrap a `pandas
<http://pandas.pydata.org/>`_ DataFrame instance, which is accessible through
the ``.data`` attribute and can be used for any manipulations that aren't
already provided by methods in the wrapper class.

``cnary``
~~~~~~~~~

.. automodule:: cnvlib.cnary
    :members:
    :undoc-members:
    :show-inheritance:

``vary``
~~~~~~~~

.. automodule:: cnvlib.vary
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_commands:

Interface to CNVkit sub-commands
--------------------------------

``commands``
~~~~~~~~~~~~

The public API for each of the commands defined in the CNVkit workflow.

.. automodule:: cnvlib.commands
    :members:
    :undoc-members:
    :show-inheritance:


The following modules implement lower-level functionality specific to each of
the CNVkit sub-commands.

``access``
~~~~~~~~~~

.. automodule:: cnvlib.access
    :members:
    :undoc-members:
    :show-inheritance:

``antitarget``
~~~~~~~~~~~~~~

.. automodule:: cnvlib.antitarget
    :members:
    :undoc-members:
    :show-inheritance:

``autobin``
~~~~~~~~~~~

.. automodule:: cnvlib.autobin
    :members:
    :undoc-members:
    :show-inheritance:

``batch``
~~~~~~~~~

.. automodule:: cnvlib.batch
    :members:
    :undoc-members:
    :show-inheritance:


``call``
~~~~~~~~

.. automodule:: cnvlib.call
    :members:
    :undoc-members:
    :show-inheritance:

``coverage``
~~~~~~~~~~~~

.. automodule:: cnvlib.coverage
    :members:
    :undoc-members:
    :show-inheritance:

``diagram``
~~~~~~~~~~~

.. automodule:: cnvlib.diagram
    :members:
    :undoc-members:
    :show-inheritance:

``export``
~~~~~~~~~~

.. automodule:: cnvlib.export
    :members:
    :undoc-members:
    :show-inheritance:

``fix``
~~~~~~~

.. automodule:: cnvlib.fix
    :members:
    :undoc-members:
    :show-inheritance:

``heatmap``
~~~~~~~~~~~

.. automodule:: cnvlib.heatmap
    :members:
    :undoc-members:
    :show-inheritance:

``importers``
~~~~~~~~~~~~~

.. automodule:: cnvlib.importers
    :members:
    :undoc-members:
    :show-inheritance:

``import_rna``
~~~~~~~~~~~~~~

.. automodule:: cnvlib.import_rna
    :members:
    :undoc-members:
    :show-inheritance:


``metrics``
~~~~~~~~~~~

.. automodule:: cnvlib.metrics
    :members:
    :undoc-members:
    :show-inheritance:

``reference``
~~~~~~~~~~~~~

.. automodule:: cnvlib.reference
    :members:
    :undoc-members:
    :show-inheritance:

``reports``
~~~~~~~~~~~

Supports the sub-commands :ref:`breaks` and :ref:`genemetrics`.

.. automodule:: cnvlib.reports
    :members:
    :undoc-members:
    :show-inheritance:

``scatter``
~~~~~~~~~~~

.. automodule:: cnvlib.scatter
    :members:
    :undoc-members:
    :show-inheritance:

``segmentation``
~~~~~~~~~~~~~~~~

.. automodule:: cnvlib.segmentation
    :members:
    :undoc-members:
    :show-inheritance:

``segmetrics``
~~~~~~~~~~~~~~~~

.. automodule:: cnvlib.segmetrics
    :members:
    :undoc-members:
    :show-inheritance:

``target``
~~~~~~~~~~

.. automodule:: cnvlib.target
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_helpers:

Helper modules
--------------

``cmdutil``
~~~~~~~~~~~

.. automodule:: cnvlib.cmdutil
    :members:
    :undoc-members:
    :show-inheritance:

``core``
~~~~~~~~

.. automodule:: cnvlib.core
    :members:
    :undoc-members:
    :show-inheritance:

``descriptives``
~~~~~~~~~~~~~~~~

.. automodule:: cnvlib.descriptives
    :members:
    :undoc-members:
    :show-inheritance:

``parallel``
~~~~~~~~~~~~

.. automodule:: cnvlib.parallel
    :members:
    :undoc-members:
    :show-inheritance:

``params``
~~~~~~~~~~

Defines several constants used in the pipeline.

.. automodule:: cnvlib.params
    :members:
    :undoc-members:
    :show-inheritance:

``plots``
~~~~~~~~~

.. automodule:: cnvlib.plots
    :members:
    :undoc-members:
    :show-inheritance:

``rna``
~~~~~~~~~

.. automodule:: cnvlib.rna
    :members:
    :undoc-members:
    :show-inheritance:

``samutil``
~~~~~~~~~~~

.. automodule:: cnvlib.samutil
    :members:
    :undoc-members:
    :show-inheritance:

``segfilters``
~~~~~~~~~~~~~~

.. automodule:: cnvlib.segfilters
    :members:
    :undoc-members:
    :show-inheritance:

``smoothing``
~~~~~~~~~~~~~

.. automodule:: cnvlib.smoothing
    :members:
    :undoc-members:
    :show-inheritance:

