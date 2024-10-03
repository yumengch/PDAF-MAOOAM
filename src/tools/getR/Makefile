MODEL_DIR=model
NETCDF_INC=-I/usr/include/openmpi-x86_64 -I/usr/lib64/gfortran/modules/openmpi
NETCDF_LIB=-L/usr/lib64/openmpi/lib -lnetcdff

# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = -fopenmp $(NETCDF_INC) -O2 -cpp -I$(MODEL_DIR)
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

PDAF_OBJECTS = mod_kind_pdaf.o mod_romb_pdaf.o getR.o

all: $(PROGRAMS)

$(PROGRAMS) : $(MODEL_OBJECTS) $(PDAF_OBJECTS)
	$(FC) -o $@ $(COMPILER_FLAGS) $^ $(LD_FLAGS)

.f90.o : 
	$(FC) $(COMPILER_FLAGS) -o $*.o -c $*.f90

clean:
	$(RM) *.o *.mod $(MODEL_DIR)/*.o $(PROGRAMS)

.PHONY: clean all
