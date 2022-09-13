INCS=-I/usr/include -I/home/yumengch/NCEO/MAOOAM/PDAF_V2.0/include
LIBNC=-L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm
LIBPDAF=-L/home/yumengch/NCEO/MAOOAM/PDAF_V2.0/lib -lpdaf-d
# Define relevant compiler options for debug case and normal case for each compiler
COMPILER_FLAGS = $(INCS) -O3 -cpp # -fcheck=all -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp
LD_FLAGS = $(LIBNC) $(LIBPDAF)  -L/usr/lib  -llapack -lblas

FC=mpif90

SHELL = /bin/sh
DEBUG = false

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam

MODEL_OBJECTS = model/params.o \
				model/util.o \
				model/inprod_analytic.o \
				model/tensor_def.o \
				model/aotensor_def.o \
				model/tl_ad_tensor.o \
				model/model_def.o \
				model/integrator_def.o \
				model/rk2_integrator.o \
				model/rk4_integrator.o \
				model/rk2_tl_integrator.o \
				model/rk2_ad_integrator.o \
				model/rk4_tl_integrator.o \
				model/rk4_ad_integrator.o \
				model/stat.o

PDAF_OBJECTS = mod_kind_pdaf.o \
          mod_romb_pdaf.o \
		  mod_parallel_pdaf.o \
		  mod_filteroptions_pdaf.o \
		  mod_inflation_pdaf.o \
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
