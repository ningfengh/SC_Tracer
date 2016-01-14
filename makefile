IDIR = ./include
CC=icpc
LFLAGS=-I$(IDIR) -fopenmp
CFLAGS=-I$(IDIR) -std=c++11 -fast -fopenmp
SDIR=./source
ODIR=./obj

default: SC_tracer


SC_tracer: $(ODIR)/SC_tracer.o $(ODIR)/raytrace.o $(ODIR)/photon_map.o $(ODIR)/triangle.o $(ODIR)/file_io.o
	$(CC) $(LFLAGS) -o SC_tracer $(ODIR)/SC_tracer.o $(ODIR)/raytrace.o $(ODIR)/photon_map.o $(ODIR)/triangle.o $(ODIR)/file_io.o

$(ODIR)/SC_tracer.o: $(SDIR)/SC_tracer.cpp
	$(CC) $(CFLAGS) -c $(SDIR)/SC_tracer.cpp  -o $(ODIR)/SC_tracer.o
	
$(ODIR)/raytrace.o: $(SDIR)/raytrace.cpp
	$(CC) $(CFLAGS) -c $(SDIR)/raytrace.cpp -o $(ODIR)/raytrace.o

$(ODIR)/photon_map.o: $(SDIR)/photon_map.cpp
	$(CC) $(CFLAGS) -c $(SDIR)/photon_map.cpp -o $(ODIR)/photon_map.o	

$(ODIR)/triangle.o: $(SDIR)/triangle.cpp
	$(CC) $(CFLAGS) -c $(SDIR)/triangle.cpp -o $(ODIR)/triangle.o
	
$(ODIR)/file_io.o: $(SDIR)/file_io.cpp
	$(CC) $(CFLAGS) -c $(SDIR)/file_io.cpp -o $(ODIR)/file_io.o		 

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o SC_tracer
