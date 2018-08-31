objects = kdtree.o projection.o
LINKS = -lm 
run:$(objects)
	mpicc -o run $(objects) $(LINKS)
projection.o:projection.c
	mpicc -c projection.c
kdtree.o:kdtree.c
	mpicc -c kdtree.c
clean:
	rm run $(objects)