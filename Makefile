CC = gfortran
CFLAGS = -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -O3
# Uncomment the appropriate CERNLPATH variable corresponding to your install
#CERNLPATH=/usr/lib/cernlib/2006-g77/lib
#CERNLPATH=/usr/lib/x86_64-linux-gnu
#CERNLPATH=/usr/lib64/

#mathlib = -lm
#cernlib = $(CERNLPATH)/libpacklib.a $(CERNLPATH)/libkernlib.a $(CERNLPATH)/libmathlib.a

# ===========
exact_test: clean exact_sim
	~/exact_sim.exe -m 300.0E-9 -t 0.01
	~/exact_sim.exe -m 300.0E-9 -t 0.001
	~/exact_sim.exe -m 3.0E3 -t 1.0E-9

exact_d2o_sim: exact_d2o_simulation.o
	$(CC) $(CFLAGS) $^ exact_banfor_module.o uuid_module.o f90getopt.o get_parameters_module.o -o ~/exact_d2o_sim.exe 

exact_sim: exact_sim.o
	$(CC) $(CFLAGS) $^ exact_banfor_module.o uuid_module.o f90getopt.o get_parameters_module.o material_list.o -o ~/exact_sim.exe

exact_d2o_simulation.o: exact_banfor_module.mod
	$(CC) $(CFLAGS) -c d2o_simulation.f95 
	mv d2o_simulation.o $@

exact_sim.o: exact_banfor_module.mod
	$(CC) $(CFLAGS) -c main.f95
	mv main.o $@

exact_banfor_module.mod: exact_banfor_module.f95 uuid_module.f90 f90getopt.F90 get_parameters_module.f95
	$(CC) $(CFLAGS) -c uuid_module.f90 
	$(CC) $(CFLAGS) -c f90getopt.F90
	$(CC) $(CFLAGS) -c material_list.f95
	$(CC) $(CFLAGS) -c get_parameters_module.f95
	$(CC) $(CFLAGS) -c exact_banfor_module.f95 


# ===========

clean:
	rm -f *.o *.txt *.dat *.mod *.exe ~/*sim.exe
