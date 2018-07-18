FC 	:= gfortran
FFLAGS	:= -Wall -g -O2 -fopenmp -Wno-unused -fcheck=all -ffpe-trap=invalid,zero,overflow

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=array-temps,bounds,do,mem,pointer\
	-g -ffpe-trap=invalid,zero,overflow
endif

OBJS	:= m_fluid_1d.o m_config.o m_generic.o m_transport_data.o

INCDIRS	:= particle_core
LIBDIRS := particle_core
LIBS	:= particle_core

VPATH := config_fortran

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: ../streamer_1d

clean:
	$(RM) *.o *.mod ../streamer_1d

# Dependency information
m_fluid_1d.o:	m_generic.o m_transport_data.o
m_particle_1d.o:	m_generic.o
m_generic.o: m_units_constants.o m_config.o
../streamer_1d:		$(OBJS)
