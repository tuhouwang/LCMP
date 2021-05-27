.PHONY: clean all

FC = gfortran

OBJS = function_MultiLC.o MultiLC.o complex_bessel.o

target = MultiLC.exe

$(target): $(OBJS)
	$(FC) $(OBJS) liblapack.a librefblas.a -o $@

%.o: %.f90
	$(FC) -c $< -o $@
	
complex_bessel.o: complex_bessel.f
	$(FC) -c $< -o $@
    
clean:
	-rm -rf *.o *.mod *.exe tl.bin a.out
