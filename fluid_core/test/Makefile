FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O3
TESTS	:= test_m_transport_data test_m_transport_schemes test_m_time_steppers

# Need these libraries for the tests
BASEDIR	:= ../../
LIBS	:= fluid_core fosito
INCDIRS := $(addprefix $(BASEDIR), $(LIBS))
LIBDIRS := ../ $(INCDIRS)

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: 	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) *.o *.mod $(TESTS)

# Dependency information
