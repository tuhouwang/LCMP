.PHONY: clean all

all:
	(cd lapack-3.10.0;	make)
	(cd srcMultiLC;     make)
    
clean:
	-rm -rf *.exe tl.bin

cleanall:
	(cd lapack-3.10.0;	make clean)
	(cd srcMultiLC;	    make clean)  
