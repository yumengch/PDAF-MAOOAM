
ifndef COMPILER
  COMPILER=gfortran
endif
INCS=-I/usr/include -I/home/yumengch/NCEO/MAOOAM/PDAF_V2.0/include
LIBNC=-L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm
LIBPDAF=-L/home/yumengch/NCEO/MAOOAM/PDAF_V2.0/lib -lpdaf-d
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2 -fpp
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -O2 -Wall -cpp $(INCS)
GCC_LD_FLAGS = $(LIBNC) $(LIBPDAF)  -L/usr/lib -llapack -lblas -lm
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl -fpp
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = $(INCS) -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp
GCC_LD_FLAGS_DBG = $(LIBNC) $(LIBPDAF)  -L/usr/lib  -llapack -lblas

FC=mpif90

SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif


ifeq ($(COMPILER),gfortran)
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
else ifeq ($(COMPILER),ifort)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam

MODULES = params.mod \
          util.mod \
          inprod_analytic.mod \
          tensor_def.mod \
          aotensor_def.mod \
          tl_ad_tensor.mod \
          model_def.mod \
          integrator_def.mod \
          rk2_integrator.mod \
          rk4_integrator.mod \
          rk2_tl_integrator.mod \
          rk2_ad_integrator.mod \
          rk4_tl_integrator.mod \
          rk4_ad_integrator.mod \
          stat.mod \
          mod_kind_pdaf.mod \
          mod_parallel_pdaf.mod \
          mod_filteroptions_pdaf.mod \
          mod_inflation_pdaf.mod \
          mod_modelwriter_pdaf.mod \
          mod_model_pdaf.mod \
          mod_localization_pdaf.mod \
          mod_config_pdaf.mod \
          mod_U_pdaf.mod \
          mod_observations_pdaf.mod \
          mod_obswriter_pdaf.mod \
          mod_U_PDAFomi_pdaf.mod \
          mod_init_pdaf.mod \
          mod_assimilate_pdaf.mod


MODULE_OBJECTS = $(MODULES:.mod=.o)

all: $(PROGRAMS) 

debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all


%.mod: %.f90
	$(FC) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%.o: %.f90 $(MODULES)
	$(FC) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%: %.o $(MODULE_OBJECTS)
	$(FC) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS) 

inprod_analytic.mod: params.mod util.mod
aotensor_def.mod: tensor_def.mod inprod_analytic.mod params.mod
integrator_def.mod: model_def.mod
rk2_integrator.mod: integrator_def.mod
tl_ad_tensor.mod: aotensor_def.mod tensor_def.mod
rk2_tl_integrator.mod: rk2_integrator.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
