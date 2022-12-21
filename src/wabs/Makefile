wabs_src = wabs.f90
wabs_obj = $(patsubst %.f90, build/$(D)/%.o, $(wabs_src))
wabs_dep = $(patsubst %.f90, build/$(D)/%.d, $(wabs_src))
wabs := build/$(D)/nnp-wabs

# === Manual commands ==========================================================
CF = gfortran
DFLAGS = -O -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace -fimplicit-none
libpath= ../lib
libobj = $(libpath)/$(wildcard $(libpath)/*.o)
libmod = $(libpath)/$(wildcard $(libpath)/*.mod)

ug_obj := $(patsubst %.f90, %.o, $(wabs_src))
ug_obj_parallel := $(patsubst %.f90, %.o, $(wabs_src))

default: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -g -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -g -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o ~/nnp-wabs

#parallel: $(ug_obj_parallel)
#
#$(ug_obj_parallel): $(patsubst %.o, %.f90, $@)
#	$(CF) $(DFLAGS) -fopenmp -I$(libpath) -c $(patsubst %.o, %.f90, $@)
#	$(CF) $(DFLAGS) -fopenmp -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o ~/nnp-wabs

clean:
	rm ~/nnp-wabs
	rm *.o
# === Manual commands ==========================================================

wabs: $(wabs)
$(wabs): $(wabs_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(wabs_dep)
