.. SpM documentation master file, created by
   sphinx-quickstart on Thu Aug 10 10:08:31 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calculation flow
===============================
In this chapter, the calculation flow in SpM is shown.


Overview of calculation flow
----------------------------
The main function for SpM is defined in G2spectrum.cpp
and the calculation is done by following procedures.

.. blockdiag::
   
   blockdiag {
       "Set initial condition" -> "Solve equation" -> "Output results";
       group{
         color ="orange"
         "Set initial condition"
       }
       group{
         color ="red"
         "Solve equation"
       }
       group{
         color ="blue"
         "Output results"
       }
   }

1. Set initial condition

   .. blockdiag::
      
      blockdiag {
          "Read parameters" -> "Read Gtau" -> "Make a kernel";

          group{
            color ="orange"
            "Read parameters" -> "Read Gtau" -> "Make a kernel";
          }
      }
   
   - Read parameters (ReadParam function in set_initial.cpp).

   - Read Gtau (read_Gtau function in G2spectrum.cpp).
     
   - Make a kernel (mesh_linear, mesh_log and MakeKernelLinear functions in kernel.cpp).
   
2. Solve equation


   .. blockdiag::
      
      blockdiag {
          "Set parameters" -> "Solve equation by \n ADMM method" -> "Get results ";

          group{
            color ="red"
            "Set parameters" -> "Solve equation by \n ADMM method" -> "Get results ";
          }
      }
   

   
   - Set parameters (SetParameters and SetFlags functions in spm_core.cpp).
     
   - Solve equation by ADMM method (SolveEquation function in spm_core.cpp).
  
     
   - Get results (GetSpectrum and GetResults functions in spm_core.cpp)

3. Output results

   This procedure is directly implemented in G2spectrum.cpp.

Sparse modeling procedure
-------------------------
to be updated...
