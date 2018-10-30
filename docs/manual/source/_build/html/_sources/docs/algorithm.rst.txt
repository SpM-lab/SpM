.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Algorithm
===============================

``SpM`` program solves the linear equation :math:`\bm{G}=K\bm{\rho}` with respect to :math:`\bm{\rho}` for given :math:`\bm{G}`.
Because of **ill-conditioned** nature of the matrix :math:`K`, a simple treatment of this equation is numerically unstable.
For example, the solution using the `Moore-Penrose pseudo-inverse matrix
<https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse>`_
results in NaN.
Even if one manages to derive a definite solution, it is quite sensitive to numerical noise and often breaks preconditions that any physical spectra must satisfy.
This becomes particularly problematic when :math:`\bm{G}` is evaluated by quantum Monte Carlo technique.

``SpM`` provides a **physical** solution which fulfills the equation of concern within a certain accuracy.
The solution satisfies the constraints such as sum rule and nonnegativity.
The engine of ``SpM`` program uses the method of **L1-norm regularization** to separate relevant information in :math:`\bm{G}` from irrelevant one which makes the spectrum unphysical. This process is automatically done without hand-tuning parameters.

For details, see the original article

    J. Otsuki, M. Ohzeki, H. Shinaoka, K. Yoshimi,
    "Sparse modeling approach to analytical continuation of imaginary-time quantum Monte Carlo data"
    `Phys. Rev. E 95, 061302(R) (2017). <https://doi.org/10.1103/PhysRevE.95.061302>`_
