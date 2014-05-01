FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O3
OBJS	:= m_fluid_dd_1d.o m_model_choice.o m_output_1d.o m_init_cond_1d.o\
	m_efield_1d.o m_particle_1d.o

INCDIRS	:= ../../fosito ../../particle_core ../../fluid_core
LIBDIRS := ../../fosito ../../particle_core ../../fluid_core
LIBS	:= part_core fluid_core fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: ../streamer_1d

clean:
	$(RM) -f *.o *.mod ../streamer_1d

# Dependency information
m_fluid_dd_1d.o: 	m_efield_1d.o m_init_cond_1d.o m_model_choice.o
m_output_1d.o: 		m_model_choice.o m_fluid_dd_1d.o m_particle_1d.o
m_particle_1d.o: 	m_init_cond_1d.o m_efield_1d.o
../streamer_1d:		$(OBJS)
