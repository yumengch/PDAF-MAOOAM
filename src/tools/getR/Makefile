MODEL_DIR=../../fortran/model
NETCDF_INC=-I/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/include
NETCDF_LIB=-L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-c-4.8.1-okwvugrilg2h3nnql5i6bhlz6mz24swi/lib -lnetcdf -L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/lib -lnetcdff -ldl -lm

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) -O2 -cpp -I../../fortran -I../../fortran/model
LD_FLAGS = $(NETCDF_LIB)

FC=mpif90

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod

PROGRAMS = getR 

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

PDAF_OBJECTS = ../../fortran/mod_kind_pdaf.o ../../fortran/mod_romb_pdaf.o getR.o

all: $(PROGRAMS)

$(PROGRAMS) : $(MODEL_OBJECTS) $(PDAF_OBJECTS)
	$(FC) -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod $(MODEL_DIR)/*.o $(PROGRAMS)

.PHONY: clean all
