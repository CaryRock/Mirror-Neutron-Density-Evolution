mult_step_src = mult_step.f90
mult_step_obj = $(patsubst %.f90, build/$(D)/%.o, $(mult_step_src))
mult_step_dep = $(patsubst %.f90, build/$(D)/%.d, $(mult_step_src))
step := build/$(D)/nnp-step

# === Manual commands ==========================================================
CF = gfortran
DFLAGS = -march=x86-64 -O3 -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace -fimplicit-none -Wno-unused-variable -Wno-compare-reals -Wno-unused-function -fdefault-real-8

DPARALLEL := -fopenmp $(DFLAGS)

libpath= ../../lib
libobj = $(wildcard $(libpath)/*.o)
libmod = $(wildcard $(libpath)/*.mod)

#ug_obj := $(patsubst %.f90, %.o, $(mult_step_src))
ug_obj_par := $(patsubst %.f90, %.o, $(mult_step_src))

#serial: $(ug_obj)

parallel: $(ug_obj_par)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o ~/nnp-mult_step_v1

$(ug_obj_par): $(patsubst &.o, %.f90, $@)
	$(CF) $(DPARALLEL) -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DPARALLEL) -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o ~/nnp-mult_step_v1

clean:
	rm *.o
	rm ~/nnp-mult_step_v1
# === Manual commands ==========================================================

step: $(step)
$(step): $(mult_step_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(mult_step_dep)
