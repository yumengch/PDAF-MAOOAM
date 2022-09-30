INCS=-I/usr/include -I$(HOME)/NCEO/MAOOAM/PDAF_V2.0/include # -I/home/users/ia923171/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/netcdf-fortran-4.5.4-arvaui536gqcge2mjsqpngr4dxzxeyxj/include
LIBNC=-L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm # -L/home/users/ia923171/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/netcdf-fortran-4.5.4-arvaui536gqcge2mjsqpngr4dxzxeyxj/lib -lnetcdff -lnetcdf -lnetcdf
LIBPDAF=-L$(HOME)/NCEO/MAOOAM/PDAF_V2.0/lib -lpdaf-d
# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(INCS) -O0 -g -cpp -fcheck=all -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp
LD_FLAGS = $(LIBNC) $(LIBPDAF)  -L/usr/lib  -llapack -lblas

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
