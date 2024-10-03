MODEL_DIR=model

PDAF_INC=-I/home/users/ia923171/PDAF_V2.2.1/include 
PDAF_LIB=-L/home/users/ia923171/PDAF_V2.2.1/lib/ -lpdaf-var 

MKLROOT = /home/users/ia923171/.conda/envs/maooam/lib
LAPACK_LIB=-Wl,--start-group $(MKLROOT)/libmkl_intel_lp64.a $(MKLROOT)/libmkl_sequential.a $(MKLROOT)/libmkl_core.a -Wl,--end-group -lm -ldl -lpthread

NETCDF_INC=-I/usr/include/mpich-x86_64 -I/usr/lib64/gfortran/modules/mpich
NETCDF_LIB=-L/usr/lib64/mpich/lib -lnetcdff

MPI_INC = -I/usr/lib64/gfortran/modules/mpich
MPI_LIB = -L/usr/lib64/mpich/lib -lmpifort -lmpi

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) $(MPI_INC) -cpp -O2 -I$(MODEL_DIR) -fPIC
# COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) $(MPI_INC) -I$(MODEL_DIR) -Wall -cpp -O0 -g -fcheck='all' -fbounds-check -fbacktrace
LD_FLAGS = $(PDAF_LIB) $(LAPACK_LIB) $(NETCDF_LIB) $(MPI_LIB)

F2PY=/home/users/ia923171/.conda/envs/maooam/bin/f2py
FC=gfortran
LD=$(FC)

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = libmaooam.so

MODEL_OBJECTS = mod_romb_pdaf.o \
                $(MODEL_DIR)/params.o \
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
				$(MODEL_DIR)/stat.o \
				mod_maooam.o

all: $(PROGRAMS)

libmaooam.so : $(MODEL_OBJECTS)
	$(FC) -shared -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)
	$(F2PY) -c mod_model.f90 -I/home/users/ia923171/PDAF-MAOOAM/src/model -L/home/users/ia923171/PDAF-MAOOAM/src/model -lmaooam -m mod_model

.F90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.F90

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod model/*.o $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out *.cpython*.so

.PHONY: clean all test %.test
