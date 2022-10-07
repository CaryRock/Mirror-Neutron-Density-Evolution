include build_rules.mk

# Make sure we can access the various library include files directly, if needed
CFLAGS += -Ilib

default: genMthetas

#all: lib coord velmesh genMthetas

#lib: lib

#coord: coord

#velmesh: velmesh

genMthetas: genMthetas

SUBDIRS=src/deltaMthetas 
#src/lib src/coord src/velmesh src/deltaMthetas

define INCLUDE_FILE
	D = $S
	include $S/Makefile
endef

$(foreach S,$(SUBDIRS),$(eval $(INCLUDE_FILE)))

.PHONY: clean
clean: $(CLEAN)
	echo "Cleaning build/..."
	rm -rf build/
	-rmdir build 2>/dev/null
