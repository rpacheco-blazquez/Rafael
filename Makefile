FCFLAGS = -O2 -fopenmp
FC      = gfortran 

ALL = test_hw_openmp  

all: $(ALL)

test_hw_openmp: test_hw_openmp.o
	$(FC) $(FCFLAGS) $^ -o $@

%.o:%.f90
	$(FC) $(FCFLAGS) -c $^

.PHONY: clean
clean:
	rm -f *.o *~
	rm -f $(ALL)
