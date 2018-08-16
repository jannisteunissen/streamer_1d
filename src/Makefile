# Compilation settings
FC 	:= gfortran
FFLAGS	:= -Wall -g -O2 -fopenmp

# Type make DEBUG=1 to do a debug build
ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=array-temps,bounds,do,mem,pointer\
	-g -ffpe-trap=invalid,zero,overflow
endif

OBJS := m_fluid_1d.o m_particle_1d.o m_config.o m_generic.o m_transport_data.o	\
m_init_cond.o m_flux_scheme.o

INCDIRS	:= particle_core
LIBDIRS := particle_core
LIBS := particle_core
VPATH := config_fortran

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: ../streamer_1d

clean:
	$(RM) *.o *.mod ../streamer_1d

# How to compile the particle_core library
particle_core/libparticle_core.a:
	$(MAKE) -C particle_core

# Recompile libparticle_core.a when source files change
particle_core/libparticle_core.a: $(wildcard particle_core/*.f90)

# Dependency information
$(OBJS): particle_core/libparticle_core.a
m_fluid_1d.o: m_generic.o m_transport_data.o
m_particle_1d.o: m_generic.o
m_generic.o: m_units_constants.o m_config.o m_init_cond.o m_flux_scheme.o
m_init_cond.o: m_config.o
streamer_1d.o: m_generic.o m_fluid_1d.o m_particle_1d.o
../streamer_1d: $(OBJS)
