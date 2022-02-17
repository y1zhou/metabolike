Installation
=====================

Prerequisites
-------------

You will need a running Neo4j database and login credentials of a user with write access to the database.
To populate the database, you will also need to acquire the `BioCyc data files`_ and the `BRENDA text file`_.
The ``metabolike`` package doesn't include the data sources since they are large and require license agreements. Please downloaded the data files and extract them for later use.

.. _`BioCyc data files`: https://biocyc.org/download.shtml
.. _`BRENDA text file`: https://www.brenda-enzymes.org/download_brenda_without_registration.php

The MetaCyc (BioCyc) database is imported to the graph database using the provided SBML file and various ``.dat`` files.
The BRENDA text file is then mapped onto the graph using common EC numbers and KEGG reaction IDs.

Latest stable version
---------------------

To install the latest stable version of the package, run the following command::

    pip install metabolike


Development version
-------------------

As the package is still in its early stages, it is recommended to get the latest version from GitHub::

    git clone --depth 1 https://github.com/y1zhou/metabolike

To install the package, run the following command:

.. code-block:: bash

    $ cd metabolike
    $ pip install .

You can also install the package using flit_::

    $ flit install

.. _flit: https://flit.readthedocs.io/en/latest/
