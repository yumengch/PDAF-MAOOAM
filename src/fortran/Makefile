MODEL_DIR=

PDAF_INC=-I/usr/include -I$(HOME)/NCEO/MAOOAM/PDAF_V2.0/include 
PDAF_LIB=-L$(HOME)/NCEO/MAOOAM/PDAF_V2.0/lib -lpdaf-d
NETCDF_INC=-I/home/users/ia923171/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/netcdf-fortran-4.5.4-arvaui536gqcge2mjsqpngr4dxzxeyxj/include
NETCDF_LIB=-L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm 


# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) -O0 -g -cpp -fcheck=all -Wall -Wextra -Wconversion -pedantic -Wno-unused-parameter -ffpe-trap=zero,overflow,underflow -cpp
LD_FLAGS = $(NETCDF_LIB) $(PDAF_LIB)  -L/usr/lib  -llapack -lblas

FC=mpif90
LD=$(FC)

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam

MODEL_OBJECTS = $(MODEL_DIR)/params.o \
				$(MODEL_DIR)/util.o \
				$(MODEL_DIR)/inprod_analytic.o \
				$(MODEL_DIR)/tensor_def.o \
				$(MODEL_DIR)/aotensor_def.o \
				$(MODEL_DIR)/tl_ad_tensor.o \
				$(MODEL_DIR)/model_def.o \
				$(MODEL_DIR)/integrator_def.o \
				$(MODEL_DIR)/rk2_integrator.o \
				$(MODEL_DIR)/rk4_integrator.o \
				$(MODEL_DIR)/rk2_tl_integrator.o \
				$(MODEL_DIR)/rk2_ad_integrator.o \
				$(MODEL_DIR)/rk4_tl_integrator.o \
				$(MODEL_DIR)/rk4_ad_integrator.o \
				$(MODEL_DIR)/stat.o

PDAF_OBJECTS = mod_kind_pdaf.o \
          mod_romb_pdaf.o \
		  mod_parallel_pdaf.o \
		  mod_filteroptions_pdaf.o \
		  mod_inflation_pdaf.o \
		  mod_statewriter_pdaf.o \
		  mod_modelwriter_pdaf.o \
		  mod_model_pdaf.o \
		  mod_localization_pdaf.o \
		  mod_config_pdaf.o \
		  mod_U_pdaf.o \
		  mod_observations_pdaf.o \
		  mod_obswriter_pdaf.o \
		  mod_U_PDAFomi_pdaf.o \
		  mod_init_pdaf.o \
		  mod_assimilate_pdaf.o \
		  maooam.o

all: $(PROGRAMS)

maooam : $(MODEL_OBJECTS) $(PDAF_OBJECTS)
	$(FC) -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)

.F90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.F90

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod model/*.o $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
