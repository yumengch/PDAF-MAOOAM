MODEL_DIR=model

PDAF_INC=-I/home/users/ia923171/PDAF_V2.2.1/include 
PDAF_LIB=-L/home/users/ia923171/PDAF_V2.2.1/lib/ -lpdaf-var 

MKLROOT = /home/users/ia923171/.conda/envs/maooam/lib
LAPACK_LIB=-Wl,--start-group $(MKLROOT)/libmkl_intel_lp64.a $(MKLROOT)/libmkl_sequential.a $(MKLROOT)/libmkl_core.a -Wl,--end-group -lm -ldl -lpthread

NETCDF_INC=-I/usr/include/mpich-x86_64 -I/usr/lib64/gfortran/modules/mpich
NETCDF_LIB=-L/usr/lib64/mpich/lib -lnetcdff

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) -cpp -O2 -I$(MODEL_DIR)
# COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) -Wall -cpp -O0 -g -fcheck='all' -fbacktrace
LD_FLAGS = $(PDAF_LIB) $(LAPACK_LIB) $(NETCDF_LIB)

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
          mod_nfcheck_pdaf.o \
		  mod_filteroptions_pdaf.o \
		  mod_inflation_pdaf.o \
		  mod_statewriter_pdaf.o \
		  mod_modelwriter_pdaf.o \
		  mod_model_pdaf.o \
		  mod_statevector_pdaf.o \
		  mod_config_pdaf.o \
          mod_U_pdaf.o \
		  mod_init_pdaf.o \
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
