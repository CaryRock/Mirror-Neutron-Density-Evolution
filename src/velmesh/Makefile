CC = gfortran
CFLAGS = -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -O3 -flto
# Uncomment the appropriate CERNLPATH variable corresponding to your install
#CERNLPATH=/usr/lib/cernlib/2006-g77/lib
#CERNLPATH=/usr/lib/x86_64-linux-gnu
#CERNLPATH=/usr/lib64/

#mathlib = -lm
#cernlib = $(CERNLPATH)/libpacklib.a $(CERNLPATH)/libkernlib.a $(CERNLPATH)/libmathlib.a

# ===========
exact_sim: exact_sim.o
	$(CC) $(CFLAGS) $^ exact_banfor_module.o file_writing.o uuid_module.o f90getopt.o get_parameters_module.o material_list.o -o ~/velmesh

exact_sim.o: exact_banfor_module.mod
	$(CC) $(CFLAGS) -c mesh_analyze.f95
	mv mesh_analyze.o $@

exact_banfor_module.mod: file_writing.f95 exact_banfor_module.f95 uuid_module.f90 f90getopt.F90 get_parameters_module.f95
	$(CC) $(CFLAGS) -c uuid_module.f90 
	$(CC) $(CFLAGS) -c f90getopt.F90
	$(CC) $(CFLAGS) -c material_list.f95
	$(CC) $(CFLAGS) -c get_parameters_module.f95
	$(CC) $(CFLAGS) -c exact_banfor_module.f95
	$(CC) $(CFLAGS) -c file_writing.f95 

# ===========

clean:
	rm -f *.o *.txt *.dat *.mod *.exe ~/*sim.exe
