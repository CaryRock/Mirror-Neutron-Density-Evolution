track_mult_step_src = track_mult_step.f90
track_mult_step_obj = $(patsubst %.f90, build/$(D)/%.o, $(track_mult_step_src))
track_mult_step_dep = $(patsubst %.f90, build/$(D)/%.d, $(track_mult_step_src))
step := build/$(D)/nnp-step

# === Manual commands =========================================================
CF = gfortran
DFLAGS := -march=x86-64 -O3 -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace
DFLAGS := $(DFLAGS) -fimplicit-none -Wno-unused-variable -Wno-compare-reals 
DFLAGS := $(DFLAGS) -Wno-unused-function -fdefault-real-8

DPARALLEL := -fopenmp $(DFLAGS)

libpath= ../lib
libobj = $(libpath)/$(wildcard $(libpath)/*.o)
libmod = $(libpath)/$(wildcard $(libpath)/*.mod)

#ug_obj := $(patsubst %.f90, %.o, $(track_mult_step_src))
ug_obj_par := $(patsubst %.f90, %.o, $(track_mult_step_src))

execDir := "$HOME/NNprime/"
name := nnp-track_mult_step

#serial: $(ug_obj)

parallel: $(ug_obj_par)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o "$(execDir)/$(name)"
#~/nnp-track_mult_step

$(ug_obj_par): $(patsubst &.o, %.f90, $@)
	$(CF) $(DPARALLEL) -I$(libpath) -c $(patsubst %.o, %.f90, $@)
	$(CF) $(DPARALLEL) -I$(libpath) $(libobj) $(patsubst %.o, %.f90, $@) -o "$(execDir)/$(name)"
#~/nnp-track_mult_step

tidy:
	rm -f *.mod;
	rm -f *.o

clean:
	make tidy
	rm -f ~/NNprime/nnp-track_mult_step
# === Manual commands ==========================================================

step: $(step)
$(step): $(track_mult_step_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(track_mult_step_dep)
