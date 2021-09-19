.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpM's documentation!
===============================
This is a documentation of Sparse Modeling (SpM) tool for analytical continuation.

What is SpM ?
-------------
A sparse-modeling tool for computing the spectral function from the imaginary-time Green function. It removes statistical errors in quantum Monte Carlo data, and performs a stable analytical continuation. The obtained spectral function fulfills the non-negativity and the sum rule. The computation is fast and free from tuning parameters.

License
--------------
This package is distributed under GNU General Public License version 3 (GPL v3).

When you publish a publication that includes results obtained using this package,
we kindly ask you to cite the following articles

- SpM method

   J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,
   "*Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data*"
   `Phys. Rev. E 95, 061302(R) (2017) <https://doi.org/10.1103/PhysRevE.95.061302>`_.

- SpM package

   K. Yoshimi, J. Otsuki, Y. Motoyama, M. Ohzeki, and H. Shinaoka,
   "*SpM: Sparse modeling tool for analytic continuation of imaginary-time Green’s function*"
   `Comput. Phys. Commun. 244, 319-323 (2019) <https://www.sciencedirect.com/science/article/abs/pii/S0010465519302103>`_.

- SpM-Pade method (if you use)

   Y. Motoyama, K. Yoshimi, and J. Otsuki,
   "*Robust analytic continuation combining the advantages of the sparse modeling approach and Pade approximation*",
   arXiv:2109.XXXXX.


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

