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

1. Create an empty directory (named ``spm.build`` in the following) and move into it:

		$ mkdir spm.build && cd spm.build

1. Call cmake:

		$ cmake ../spm.src

1. Compile the code:

		$  make

Then, the executable file ``SpM.out``  is created in the ``spm.build/src`` directory.

## Running samples

Sample data are provided both for fermionic and bosonic cases in directories ``samples/fermion/`` and ``samples/boson/``, respectively.
Here, explanation is given for the fermionic case.

User can execute all the procedure explained below using a single script file ``run.sh``.
Edit the first two variables ``file_exe="../SpM.out"`` and ``dir_plt="../plt"`` according to your system configuration, and then run the script by

    $ ./run.sh

If succeeded, text files containing numerical data and graphs in pdf format are created in directory ``output``.


### How to generate PDFs

User can generate graphs using gnuplot.
Move into directory ``output``, and type

    gnuplot path_to_SpM/samples/plt/*

Next, move into directory ``lambda_opt`` and type

    gnuplot path_to_SpM/samples/plt/lambda_fix/*

Then some graphs such as *spectrum.pdf*, *find_lambda_opt.pdf* are generated. 
See the documantation in the homepage of the SpM for more deteils.

## Official page
The official page of the SpM is [here](xxx) (not linked yet).  