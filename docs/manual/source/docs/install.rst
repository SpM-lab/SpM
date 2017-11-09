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


Sample scripts
------------------

Some sample data are provided in ``samples`` directory:

- ``samples/fermion``  # sample for fermionic spectrum (data in the article)
- ``samples/boson``  # sample for bosonic spectrum

A script file is also provided to run through the program.
Enter into one of the sample directory, and execute the script by

::

    $ ./run.sh

You may need to change the path to ``SpM.out`` in the script.
If succeeded, results including graphs in pdf format are created in ``output`` directory.
For details of the sample script, see :ref:`tutorials`.



..
  1. Parameter file

     A simple way is to give parameters in a text file, and pass it by -i option:

     ``$  spm.build/src/SpM.out -i param.in``

     Typical input is given in ``param.in`` in the sample directories.
     A default value is used for parameters not given in the file.
     The output data will be created in ``output`` directory.

  2. Command-line arguments

     Alternatively, you can pass **all** parameters to the program as command-line arguments.
     See the script file ``run.sh`` for details.
     You can run it simply by

     ``$ ./run.sh``

     Here, it is assumed that the executable ``SpM.out`` is located in ``samples`` directory.
     If not, copy or link the executable or modify ``run.sh``.
     In the script, gnuplot is called after calculations to generate pdf files.
