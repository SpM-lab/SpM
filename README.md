SPM
====

Sparse Modeling tool for analytic continuation.
The algorithm is written in the article
"Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data"
[Phys. Rev. E 95, 061302(R) (2017).](https://doi.org/10.1103/PhysRevE.95.061302)

## Licence

This package is distributed under GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).

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

Some sample data are prepared in ``samples`` directory:

* samples/fermion  # sample for fermionic spectrum (data in the article)
* samples/boson  # sample for bosonic spectrum

There are two ways to run the program.

1. Parameter file:

	A simple way is to give parameters in a text file, and pass it by -i option:

		$  spm.build/src/SpM.out -i param.in

	Typical input is given in ``param.in`` in the sample directories.
	A default value is used for parameters not given in the file.
	The output data will be created in ``output`` directory.

1. Command-line arguments

	Alternatively, you can pass **all** parameters to the program as command-line arguments. See the script file ``run.sh`` for details.
	You can run it simply by

		$ ./run.sh

	Here, it is assumed that the executable ``SpM.out`` is located in ``samples`` directory. If not, copy or link the executable or modify ``run.sh``.
	In the script, gnuplot is called after calculations to generate pdf files.

## Input data

1. Paramete file (param.in)

	- INPUT/OUTPUT
		- statistics	: [String]	Choose "fermion" or "boson" (Default: "fermion").
		- beta	: [Double]	The inverse temperature (Default: 100).
		- filein_G	: [String]	The name of an input file for Green's function (Default: "Gtau.in").
		- column	: [Integer]	The column number where the values of G(tau) are stored in the ``filein_G`` file (Default: 1).
		- fileout_spec	: [String]	The name of output file of spectrum (Default: "spectrum.out").

	- OMEGA
		- Nomega	:	[Integer]	The number of omega to be calculated (Default: 1001).
		- omegamin	:	[Double]	The minimum value of omega (Default: -4).
		- omegamax	:	[Double]	The maximum value of omega (Default: 4).

		Note: The i-th omega, Omega[i] (i=0-Nomega), is given by omegamin+(omegamax-omegamin)/NOmega * i.

	- SVD
		- SVmin	:	[Double] Truncation value of singlular values (Default: 1e-10).  

	- ADMM
		- lambdalogmesh	:	[Double]	The log mesh of lambda (Default: 0.2).
		- lambdalogbegin :	[Double]	The log value of maximum lambda, i.e. lambda_max is given by 10^{lambdalogbegin} (Default: 0)
		- lambdalogend	:	[Double]	 The log value of minimum lambda, i.e. lambda_min is given by 10^{lambdalogend}	(Default: -1)
		- Nlambda	:	[int] The number of lambda to be calculated.
		- penalty	:	[Double]  The value of penalty coefficient. If negative, penalty is optimized during the iteration starting with its absolute value (Default: 10.0).
		- tolerance	:	[Double] The criteria of convergience (Default: 1e-6).
		- maxiteration	:	[Integer]	The maximum number of iterations (Default: 1000).
		- printlevel	:	[Integer]	0; minimum, 1; moderate, 2; verbose (Default: 2).

2. Green's function (Gtau.in)  
In SPM, the values of Green's function is only used for calculation, i.e. tau is automatically determined by the beta and the step. Please indicate the column number where the values of G(tau) are stored by "column" in the prameter file.


## Output data
Calculated results are produced in directory 'output'.
The list of files and brief explanations are given below:

```
output
├── find_lambda_opt.dat  : finding the optimal value of lambda [Fig. 4(b) in PRE]
├── lambda_dep.dat       : lambda dependence of the square error, etc. [Fig. 4(a)]
├── spectrum.dat         : spectrum for optimal value of lambda [Fig. 2(b2)]
├── SV.dat               : singular values
├── SVD_U.dat            : (optional)
├── SVD_V.dat            : (optional)
├── mesh_omega.dat       : (optional)
├── mesh_tau.dat         : (optional)
├── lambda_opt
│   ├── iter.dat     : convergence of solution
│   ├── x_sv.dat     : $\rho_l$ [Fig. 3(b)]
│   ├── x_tw.dat     : $\rho(\omega)$ [Fig. 2(b)]
│   ├── y_sv.dat     : $G_l$ [Fig. 3(a)]
│   └── y_tw.dat     : $G(\tau)$ [Fig.2(a)]
├── lambda
│   ├── lambda_1.00e+00  : (optional) as in lambda_opt
│   ├── ...              :            for all values of lambda
```

## How to generate PDFs

Gnuplot script files for generating pdf are prepared in `samples/plt`.
After calculation is done, type as follows

	$ cd output
	$ gnuplot path_to_plt/*

More detailed plots for a fixed lambda can be generated by

	$ cd lambda_opt
	$ gnuplot path_to_plt/lambda_fix/*
