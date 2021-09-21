# SpM

Sparse Modeling tool for analytical continuation from the imaginary-time Green's function to real-requency spectral function.

The algorithm is presented in the article

* J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,  
"Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data"  
[Phys. Rev. E 95, 061302(R) (2017).](https://doi.org/10.1103/PhysRevE.95.061302)

## Licence

This package is distributed under GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).

We kindly ask you to cite the below articles
in publications that include results obtained using this package.

* The article about original algorithm
    * J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,  
    "Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data"  
    [Phys. Rev. E 95, 061302(R) (2017).](https://doi.org/10.1103/PhysRevE.95.061302)
* The article about this package
    * K. Yoshimi, J. Otsuki, Y. Motoyama, M. Ohzeki, and H. Shinaoka, "SpM: Sparse modeling tool for analytic continuation of imaginary-time Green's function" [Comput. Phys. Commun. 244, 319-323 (2019)](https://www.sciencedirect.com/science/article/pii/S0010465519302103).
* The article about the SpM-Pade method (please cite if you use)
    * Y. Motoyama, K. Yoshimi, and J. Otsuki, "Robust analytic continuation combining the advantages of the sparse modeling approach and Pade approximation" [arXiv:2109.08370](https://arxiv.org/abs/2109.08370)

## Authors

Junya Otsuki, Kazuyoshi Yoshimi, Yuichi Motoyama, Hiroshi Shinaoka, Masayuki Ohzeki

## Requirement

* LAPACK, BLAS
* FFTW3
* cpplapack (included in this package)

## How to build

### Getting the source codes

Download the latest source codes by

    $ git clone https://github.com/j-otsuki/SpM.git spm.src

Then, the source codes are downloaded in the directory `spm.src`.

### Using Cmake

Build with cmake command is done in a separate directory, e.g. `spm.build`.
Type the following commands:

    $ mkdir spm.build && cd spm.build
    $ cmake ../spm.src
    $ make

Then, the executable file `SpM.out` is created in directory `spm.build/src`.

## Sample scripts

Some sample data are provided in `samples` directory:

* `samples/fermion`  # sample for fermionic spectrum (data in the article)
* `samples/boson`  # sample for bosonic spectrum

A script file, `run.sh`, is provided to run through the program.
Enter into the directory `samples/fermion`, and execute the script by

    $ ./run.sh

You may need to change the path to `SpM.out` in the script.
If succeeded, results including graphs in eps format are created in `output` directory.
For details, see the document linked below.

## Official page

The official page of the SpM is [here](https://spm-lab.github.io/SpM/manual/build/html/index.html).  
