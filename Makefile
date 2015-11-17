FC 		:= gfortran
FFLAGS	:= -Wall -Wextra -O2
OBJS	:= m_transport_data.o m_time_steppers.o m_transport_schemes.o
INCDIRS :=

ifeq ($(DEBUG), 1)
	FFLAGS += -O0 -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
			-pedantic -finit-real=snan
endif

ifeq ($(PROF), 1)
	FFLAGS += -pg
endif

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

.PHONY: all test clean

all: 	libfluid_core.a

libfluid_core.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

test: 	libfluid_core.a
	$(MAKE) -C test

clean:
	$(RM) *.o *.mod
	$(MAKE) -C test clean

# Dependency information
