lib_src = material_list.f90 f90getopt.f90 file_writing.f90 get_parameters_module.f90 uuid_module.f90 distributions.f90 exact_banfor_module.f90
# These could probably benefit from a couple judicious `foreach` calls

lib_obj = $(patsubst %.f90, build/$(D)/%.o, $(lib_src))
lib_dep = $(patsubst %.f90, build/$(D)/%.o, $(lib_src))
lib	:=build/$(D)/lib.a

# === MANUAL ===================================================================
ug_obj  = $(patsubst %.f90, %.o, $(lib_src))
DFLAGS = -march=x86-64 -g -cpp -Wall -Wextra -fbacktrace -fcheck=all -pedantic -fimplicit-none -fdefault-real-8

all: $(ug_obj)

$(ug_obj): $(patsubst %.o, %.f90, $@)
	gfortran $(DFLAGS) -g $(patsubst %.o, %.f90, $@) -c

clean:
	rm *.o *.mod
# === MANUAL ===================================================================	

lib: $(lib)
$(lib): $(lib_obj)
	echo " AR $@..."
	$(AR) -cr $@ $^

-include $(lib_dep)
