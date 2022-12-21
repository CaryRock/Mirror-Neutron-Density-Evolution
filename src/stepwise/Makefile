step_src = step.f90
step_obj = $(patsubst %.f90, build/$(D)/%.o, $(step_src))
step_dep = $(patsubst %.f90, build/$(D)/%.d, $(step_src))
step := build/$(D)/nnp-step

# === Manual commands ==========================================================
CF = gfortran
DFLAGS = -Og -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace -fimplicit-none -fdefault-real-16 -fopenmp
libpath= ../lib
libobj = $(libpath)/$(wildcard $(libpath)/*.o)
libmod = $(libpath)/$(wildcard $(libpath)/*.mod)

ug_obj := $(patsubst %.f90, %.o, $(step_src))

default: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o ~/nnp-step

clean:
	rm ~/nnp-step
	rm *.o
# === Manual commands ==========================================================

step: $(step)
$(step): $(step_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(step_dep)
