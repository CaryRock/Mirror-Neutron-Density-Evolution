gen_src = genMtheta.f95
gen_obj = $(patsubst %.f95,build/$(ARCH)/$(D)/%.o,$(gen_src))
gen_dep = $(patsubst %.f95,build/$(ARCH)/$(D)/%.dep,$(gen_src))
gen := build/$(ARCH)/$(D)/nnp-genMthetas

default: generateMtheta

generateMtheta: genMtheta.f95
	gfortran -O3 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace genMtheta.f95 -o ~/genMtheta

gen: $(gen)
$(gen): $(gen_obj)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS)

-include $(gen_dep)
