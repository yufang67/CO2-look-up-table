#
#  Makefile for CO2 table:
#
#  Type "make" or "make rp" to create .o files needed in 
#  applications directories
#
#
#
FFLAGS = -g -O0
LFLAGS = -c
F77    = gfortran -c
LINK   = gfortran

#.f.o: ; $(LINK)  $(LFLAGS) $*.f90 
 
%.o : %.f90
	$(LINK) $(LFLAGS)  $<

OBJECTS = \
  Def_constants.o \
  Def_variables.o \
  deriv_disfonc.o \
  deriv_expfonc.o \
  helmholtz_deriv.o \
  helmholtz_dimless.o \
  Properties.o \
  non_linear_solvers.o \
  grid_functions.o \
  grid_construction_LL.o \
  grid_construction_LH.o \
  grid_construction_R.o \
  grid_construction_HT.o \
  grid_construction_TPL.o \
  grid_construction_TPM.o \
  grid_construction_TPH.o \
  saturation_curve.o \
  Grid.o \
  var_const.o \
  interp_functions.o \
  Interp_table.o \
  Transprop.o \
  Derivees.o \
  solver_eos.o 

CO2SOURCES = \
  Def_constants.f90 \
  Def_variables.f90 \
  deriv_disfonc.f90 \
  deriv_expfonc.f90 \
  helmholtz_deriv.f90 \
  helmholtz_dimless.f90 \
  Properties.f90 \
  non_linear_solvers.f90 \
  grid_functions.f90 \
  grid_construction_LL.f90 \
  grid_construction_LH.f90 \
  grid_construction_R.f90 \
  grid_construction_HT.f90 \
  grid_construction_TPL.f90 \
  grid_construction_TPM.f90 \
  grid_construction_TPH.f90 \
  saturation_curve.f90 \
  Grid.f90 \
  var_const.f90 \
  interp_functions.f90 \
  Interp_table.f90 \
  Transprop.f90 \
  Derivees.f90 \
  solver_eos.f90 

CO2: $(OBJECTS)

clean:
	rm -f *.o *.mod *.MOD



### DO NOT remove this line - make depends on it ###
