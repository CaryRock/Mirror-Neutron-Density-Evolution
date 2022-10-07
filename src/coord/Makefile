coord_src = nnp-coord.f90
coord_obj = $(patsubst %.f90, build/$(D)/%.o, $(coord_src))
coord_dep = $(patsubst %.f90, build/$(D)/%.d, $(coord_src))
coord :=build/$(D)/nnp-coord

# === Workaround for gfortran -MM -MG -MT bug ==================================
DFLAGS = -O3 -cpp -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace

lib_o = $(wildcard ../lib/*.o)
lib_mod = $(wildcard ../lib/*.mod)

ug_obj = $(patsubst %.f90, %.o, $(coord_src))

default: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	cp -r ../lib/$(lib_o) ../lib/$(lib_mod) .
	gfortran $(DFLAGS) $(lib_o) $(patsubst %.o, %.f90, $@) -o $(patsubst %.o, ~/nnp-coord, $@)
	rm *.o *.mod
# === Workaround for gfortran -MM -MG -MT bug ==================================

coord: $(coord)
$(coord): $(coord_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(coord_dep)
