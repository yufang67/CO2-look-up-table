# CO2 Look Up Table
This model is basically written by FORTRAN 90. It can compute the CO2 thermodynamic properties and transport porperties in liquid, gas, supercritical and liquid-gas two-phase states. The inputs are the density (specific volume) and the internal energy. This input choice is coherent with most of density-based flow solvers. More details refer to [An Accurate and Efficient Look-up Table Equation of State for Two-Phase Compressible Flow Simulations of Carbon Dioxide](https://pubs.acs.org/doi/10.1021/acs.iecr.8b00507)


---
### - ***Code organisation***
### - ***How to use***
---

## 1. Code organisation

- ***Original Span-Wagner EoS*** (thermodynamic properties):
  - helmholtz_deriv.f90
  - helmholtz_dimless.f90
  - deriv_disfonc.f90
  - deriv_expfonc.f90
  - Def_constants.f90
  - Properties.f90
  - saturation_curve.f90
  - Derivees.f90
- ***Make look-up table***:
  - Def_variables.f90
  - grid_construction_HT.f90
  - grid_construction_LH.f90
  - grid_construction_LL.f90
  - grid_construction_R.f90
  - grid_construction_TPH.f90
  - grid_construction_TPM.f90
  - grid_construction_TPL.f90
  - grid_functions.f90
  - interp_functions.f90
  - Grid.f90
  - Interp_table.f90
- ***Root-finding solver***:
  - non_linear_solvers.f90
  - axl_solvers.f90
  -  solver_eos.f90
- ***Transport properteis***:
  -  Transprop.f90

- ***Test/evaluation cases***:
  - program/
  
  ---
  ## 2. How to use
  - Compliler: gfortran, lapack library, openblas library (***Different version can lead divergence in some points during grid construction. To fix it, need adapt initial guess. You may find a clue in grid construction modules***)
  - Test cases are in ***program*** file.
  - To couple to OpenFoam (C++ based) solver:
    - compile .f90 using Makefile ("make")
    - create .a library ("ar rc libco2lib.a *.o")
    - check library ("nm libco2lib.a")
    - compile .a library with OpenFOAM solver (setup Make/option in OpenFOAM solver)
  
  -  ***Main Module***:
    - Grid.f90 : construct table
    - Interp_table.f90 : use table to compute properties
    - Properties.f90 : use original SW EoS to compute properties
    - Transprop.f90 : compute transport properties
    - Derivees.f90 : compute derivatives of thermal variables
 
