test_src = test_exact_banfor.f90
test_obj = $(patsubst %.f90, build/$(D)/%.o, $(coord_src))
test_dep = $(patsubst %.f90, build/$(D)/%.d, $(coord_src))
test :=build/$(D)/nnp-test

test: $(test)
$(test): $(test_obj) $(lib)
	echo "LD $@..."
	$(CF) -o $@ $^ $(LFLAGS) $(lib)

-include $(test_dep)
