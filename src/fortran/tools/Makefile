MODEL_DIR=model

PDAF_INC=-I/home/users/yumengch/PDAF_V2.0/include
PDAF_LIB=-L/home/users/yumengch/PDAF_V2.0/lib -lpdaf-var -L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netlib-lapack-3.10.1-2ugsub76cuzvc4ozk5oyjs27cstb23wh/lib64 -llapack -lblas
NETCDF_INC=-I/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/include
NETCDF_LIB=-L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-c-4.8.1-okwvugrilg2h3nnql5i6bhlz6mz24swi/lib -lnetcdf -L/home/users/yumengch/spack/opt/spack/linux-centos7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-mj2hjugxr6ngqu5iwejgs6jmplozc7zb/lib -lnetcdff -ldl -lm

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) $(PDAF_INC) -O2 -cpp
LD_FLAGS = $(NETCDF_LIB) $(PDAF_LIB)
#

FC=mpif90

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = genCovar

MODEL_OBJECTS = ../model/params.o \
				../model/util.o \
				../model/inprod_analytic.o \
				../model/tensor_def.o \
				../model/aotensor_def.o \
				../model/tl_ad_tensor.o \
				../model/model_def.o \
				../model/integrator_def.o \
				../model/rk2_integrator.o \
				../model/rk4_integrator.o \
				../model/rk2_tl_integrator.o \
				../model/rk2_ad_integrator.o \
				../model/rk4_tl_integrator.o \
				../model/rk4_ad_integrator.o \
				../model/stat.o

PDAF_OBJECTS = ../mod_kind_pdaf.o \
          ../mod_romb_pdaf.o \
		  genCovar.o

all: $(PROGRAMS)

genCovar : $(MODEL_OBJECTS) $(PDAF_OBJECTS)
	$(FC) -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)

.F90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.F90

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod ../*.o ../model/*.o $(PROGRAMS) tests/*.out

.PHONY: clean all test %.test
