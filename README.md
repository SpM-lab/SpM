SPM
====

Sparse Modeling tool for analytical continuation.

The algorithm is presented in the article

* J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,  
"Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data"  
[Phys. Rev. E 95, 061302(R) (2017).](https://doi.org/10.1103/PhysRevE.95.061302)

## Licence

This package is distributed under GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).

We kindly ask you to cite the article above
in publications that include results obtained using this package.

## Author
Junya Otsuki, Kazuyoshi Yoshimi, Hiroshi Shinaoka, Masayuki Ohzeki

## Requirement

* LAPACK, BLAS
* cpplapack (included in this package)


## How to build

### Getting the source codes

Download the latest source codes by

	$ git clone https://github.com/j-otsuki/SpM.git spm.src

Then, the source codes are downloaded in the directory ``spm.src``.

### Using Cmake

Build with cmake command is done in a separate directory, e.g. ``spm.build``.
Type the following commands:

	$ mkdir spm.build && cd spm.build
	$ cmake ../spm.src
	$ make

Then, the executable file ``SpM.out`` is created in directory ``spm.build/src``.

## Sample scripts

Some sample data are provided in ``samples`` directory:

* ``samples/fermion``  # sample for fermionic spectrum (data in the article)
* ``samples/boson``  # sample for bosonic spectrum

A script file is also provided to run through the program.
Enter into one of the sample directory, and execute the script by

    $ ./run.sh

You may need to change the path to ``SpM.out`` in the script.
If succeeded, results including graphs in pdf format are created in ``output`` directory.
For details, see the document linked below.


## Official page
The official page of the SpM is [here](xxx) (not linked yet).  
