.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

How to install
===============================

Requirement
--------------
* LAPACK, BLAS
* cpplapack (included in this package)


Download
---------
* **Stable version**

  You can download a stable version of the ``SpM`` program in https://github.com/SpM-lab/SpM/releases. After download the archive, decompress it by

  ::

    $ tar zxvf v1.0.tar.gz

* **Latest codes**

  You can get the latest codes from `our GitHub repository <https://github.com/SpM-lab/SpM>`_. Execute the following command

  ::

    $ git clone https://github.com/SpM-lab/SpM.git spm.src

  and then the source codes are downloaded in the directory ``spm.src``.

Build
------

1. First, create an empty directory and move into it:

  ::

    $ mkdir spm.build && cd spm.build

  We assume that ``spm.build`` is located on the same level as ``spm.src``.

2. Call cmake

  ::

    $ cmake ../spm.src

  and all necessary files including ``Makefile`` are created according to your system configuration. [any cmake options?]

3. Now, you can compile the codes by

  ::

    $ make

  The executable file ``SpM.out`` is created in the ``spm.build/src`` directory.


Test
------------------

Some sample data are provided for test calculations.
Just enter into ``samples/fermion`` directory and execute it by

::

    $ ./run.sh

You may need to change the parameter ``file_exe="../SpM.out"`` in the script (absolute or relative path to ``SpM.out``) according to your environment.
If succeeded, results including graphs in EPS format are generated in ``output`` directory.
For details of the sample script, see :ref:`tutorials`.
