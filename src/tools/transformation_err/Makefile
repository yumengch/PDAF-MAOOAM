MODEL_DIR=../../MAOOAM-1.4/fortran

NETCDF_INC=-I/opt/cray/pe/netcdf-hdf5parallel/4.7.4.7/gnu/9.1/include
NETCDF_LIB=-L/opt/cray/pe/netcdf-hdf5parallel/4.7.4.7/gnu/9.1/lib -L/opt/cray/pe/hdf5-parallel/1.12.0.7/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz 

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(NETCDF_INC) -O2 -cpp
LD_FLAGS = $(NETCDF_LIB)

FC=ftn

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod

PROGRAMS = compare 

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

PDAF_OBJECTS = mod_kind_pdaf.o mod_romb_pdaf.o compare.o

all: $(PROGRAMS)

compare : $(MODEL_OBJECTS) $(PDAF_OBJECTS)
	$(FC) -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod $(MODEL_DIR)/*.o $(PROGRAMS)

.PHONY: clean all
