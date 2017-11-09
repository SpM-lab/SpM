.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpM's documentation!
===============================
This is a documantation of Sparse Modeling (SpM) tool for analytical continuation.

What is SpM ?
-------------
A sparse-modeling tool for computing the spectral function from the imaginary-time Green function. It removes statistical errors in quantum Monte Carlo data, and performs a stable analytical continuation. The obtained spectral function fulfills the non-negativity and the sum rule. The computation is fast and free from tuning parameters.

License
--------------
This package is distributed under GNU General Public License version 3 (GPL v3).

We kindly ask you to cite the article

    J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,
    "*Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data*"
    `Phys. Rev. E 95, 061302(R) (2017) <https://doi.org/10.1103/PhysRevE.95.061302>`_.

in publications that includes results obtained using this package.


Contents
--------
.. toctree::
   :maxdepth: 2
   :numbered:
   :glob:
   
   docs/install.rst
   docs/tutorials.rst
   docs/algorithm.rst
   docs/function.rst
   docs/inputfile.rst
   docs/outputfile.rst

