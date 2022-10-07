vel_src = mesh_analyze.f90
vel_obj = $(patsubst %.f90, build/$(D)/%.o, $(vel_src))
vel_dep = $(patsubst %.f90, build/$(D)/%.d, $(vel_src))
vel := build/$(D)/nnp-velmesh

# === Workaround for gfortran -MM -MG -MT bug ==================================
DFLAGS = -O3 -cpp -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace 

lib_o = $(wildcard ../lib/*.o)
lib_mod = $(wildcard ../lib/*.mod)

ug_obj = $(patsubst %.f90, %.o, $(vel_src))

default: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	cp -r ../lib/$(lib_o) ../lib/$(lib_mod) .
	gfortran $(DFLAGS) $(lib_o) $(patsubst %.o, %.f90, $@) -o $(patsubst %.o, ~/nnp-velmesh, $@)
	rm *.o *.mod
# === Workaround for gfortran -MM -MG -MT bug ==================================

vel: $(vel)
$(vel): $(vel_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(vel_dep)
