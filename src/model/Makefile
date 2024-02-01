MODEL_DIR=model

PDAF_INC=
PDAF_LIB=
NETCDF_INC=-I/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/include
NETCDF_LIB=-L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-c-4.8.1-okwvugrilg2h3nnql5i6bhlz6mz24swi/lib -lnetcdf -L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/lib -lnetcdff -ldl -lm

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) -cpp -O2 -fPIC
LD_FLAGS = $(NETCDF_LIB)

F2PY=/home/users/yumengch/miniconda3/envs/qgs/bin/f2py
FC=mpif90
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
	$(F2PY) -c mod_model.f90 -L. -lmaooam -m mod_model

.F90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.F90

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod model/*.o $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out *.cpython*.so

.PHONY: clean all test %.test
