scikit-genome package
=====================

Module ``skgenome`` contents
----------------------------

.. automodule:: skgenome
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_tabio:

Tabular file I/O (tabio)
------------------------

``tabio``
~~~~~~~~~

.. automodule:: skgenome.tabio
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_sk_gary:

Base class: GenomicArray
------------------------

The base class of the core objects used throughout CNVkit and scikit-genome is
GenomicArray. It wraps a `pandas <http://pandas.pydata.org/>`_ DataFrame
instance, which is accessible through the ``.data`` attribute and can be used
for any manipulations that aren't already provided by methods in the wrapper
class.

``gary``
~~~~~~~~

.. automodule:: skgenome.gary
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_sk_interval:

Genomic interval arithmetic
---------------------------

``intersect``
~~~~~~~~~~~~~

.. automodule:: skgenome.intersect
    :members:
    :undoc-members:
    :show-inheritance:

``merge``
~~~~~~~~~

.. automodule:: skgenome.merge
    :members:
    :undoc-members:
    :show-inheritance:

``subdivide``
~~~~~~~~~~~~~

.. automodule:: skgenome.subdivide
    :members:
    :undoc-members:
    :show-inheritance:

``subtract``
~~~~~~~~~~~~

.. automodule:: skgenome.subtract
    :members:
    :undoc-members:
    :show-inheritance:


.. _api_sk_helpers:

Helper modules
--------------

``chromsort``
~~~~~~~~~~~~~

.. automodule:: skgenome.chromsort
    :members:
    :undoc-members:
    :show-inheritance:

``combiners``
~~~~~~~~~~~~~

.. automodule:: skgenome.combiners
    :members:
    :undoc-members:
    :show-inheritance:

``rangelabel``
~~~~~~~~~~~~~~

.. automodule:: skgenome.rangelabel
    :members:
    :undoc-members:
    :show-inheritance:
