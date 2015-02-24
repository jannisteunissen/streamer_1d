FC 	:= gfortran
FFLAGS	:= -Wall -O2 -fopenmp

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=array-temps,bounds,do,mem,pointer\
	-g -ffpe-trap=invalid,zero,overflow
endif

OBJS	:= m_fluid_dd_1d.o m_model_choice.o m_output_1d.o m_init_cond_1d.o\
	m_efield_1d.o m_particle_1d.o m_phys_domain.o



INCDIRS	:= ../fosito ../particle_core ../fluid_core
LIBDIRS := ../fosito ../particle_core ../fluid_core
LIBS	:= particle_core fluid_core fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: ../streamer_1d

clean:
	$(RM) *.o *.mod ../streamer_1d

# Dependency information
m_fluid_dd_1d.o: 	m_efield_1d.o m_init_cond_1d.o m_model_choice.o m_phys_domain.o
m_efield_1d.o:		m_phys_domain.o
m_output_1d.o: 		m_model_choice.o m_fluid_dd_1d.o m_particle_1d.o
m_particle_1d.o: 	m_init_cond_1d.o m_efield_1d.o m_phys_domain.o
m_init_cond_1d.o:	m_phys_domain.o
../streamer_1d:		$(OBJS)
