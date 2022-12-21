vel_src = mesh_analyze.f90
vel_obj = $(patsubst %.f90, build/$(D)/%.o, $(vel_src))
vel_dep = $(patsubst %.f90, build/$(D)/%.d, $(vel_src))
vel := build/$(D)/nnp-velmesh

# === Workaround for gfortran -MM -MG -MT bug ==================================
CF = gfortran
DFLAGS = -O3 -cpp -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace 
libpath= ../lib
libobj = $(libpath)/$(wildcard $(libpath)/*.o)
libmod = $(libpath)/$(wildcard $(libpath)/*.mod)

ug_obj := $(patsubst %.f90, %.o, $(vel_src))

default: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	$(CF) $(DFLAGS) -I$(libpath) -c $^
	$(CF) $(DFLAGS) -I$(libpath) $(libobj) $^ -o nnp-velmesh
# === Workaround for gfortran -MM -MG -MT bug ==================================

vel: $(vel)
$(vel): $(vel_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(vel_dep)
