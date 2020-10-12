FC 	:= gfortran
FFLAGS	:= -Wall -g -O2 -fopenmp $(USR_FLAGS)
TESTS	:= test_m_cross_sec test_openmp test_swarm ion_drift test_outside_check

# Need these libraries for the tests
LIBS	:= particle_core
INCDIRS := ..
LIBDIRS := ..

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: 	$(TESTS)

test:	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) *.o *.mod $(TESTS)

# Dependencies
$(TESTS): ../libparticle_core.a
