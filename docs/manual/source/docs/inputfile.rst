.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _inputfiles:

Input files
===============================
1. Parameter file (param.in)
  
The variables for which default value is not given are mandatory.

* INPUT/OUTPUT

  .. csv-table::
     :header-rows: 1
     :widths: 1,1,2,4

     Name, Type, Default value, Description
     statistics, String, ---, Choose "fermion" or "boson".
     beta, Double, ---, The inverse temperature
     filein_G, String, Gtau.in, The name of an input file for Green's function.
     column, Integer, 1, The column number where the values of G(tau) are stored in the filein_G file.
     fileout_spec, String, spectrum.out, The name of output file of spectrum.
       

     
* OMEGA

  .. csv-table::
     :header-rows: 1
     :widths: 1,1,2,4

     Name, Type, Default value, Description
     NOmega, Integer, 1001, The number of omega to be calculated.
     omegamin, Double, -4, The minimum value of omega.
     omegamax, Double, 4, The maximum value of omega.

  Note: The i-th omega where i=[0:Nomega) is given by
  :math:`\verb|omega|[i]=\verb|omegamin|+(\verb|omegamax|-\verb|omegamin|)/(\verb|NOmega|-1) * i`.

* SVD

  .. csv-table::
     :header-rows: 1
     :widths: 1,1,2,4

     Name, Type, Default value, Description
     SVmin, Double, 1e-10, Truncation value of singiular values.

     
* ADMM

  .. csv-table::
     :header-rows: 1
     :widths: 1,1,2,4

     Name, Type, Default value, Description
     lambdalogmesh, Double, 0.2, The log mesh of lambda.
     lambdalogbegin, Double, 0, The log value of maximum lambda. lambda_max is given by :math:`10^{\verb|lambdalogbegin|}`.
     lambdalogend, Double, -1, The log value of minimum lambda. lambda_min is given by :math:`10^{\verb|lambdalogend|}`
     Nlambda, Integer, ---,  The number of lambda to be calculated.
     penalty, Double, 10.0, "The value of penalty coefficient. If negative, penalty is optimized during the iteration starting with its absolute value."
     tolerance, Double, 1e-6, The criteria of convergence.
     maxiteration, Integer,1000,	The maximum number of iterations.
     printlevel, Integer,2,	"0; minimum, 1; moderate, 2; verbose."

2. Green's function (Gtau.in)

   In SPM, the values of Green's function is only used for calculation,
   i.e. tau is automatically determined by the beta and the step.
   Please indicate the column number where the values of G(tau) are
   stored by "column" in the parameter file.
