# Uncomment the appropriate CERNLPATH variable corresponding to your install
#CERNLPATH=/usr/lib/cernlib/2006-g77/lib
#CERNLPATH=/usr/lib/x86_64-linux-gnu
#CERNLPATH=/usr/lib64

#mathlib = -lm
#cernlib = $(CERNLPATH)/libpacklib.a $(CERNLPATH)/libkernlib.a $(CERNLPATH)/libmathlib.a

# ==========
ARCH?=$(shell uname -m)

CF	= gfortran
AR	= ar
STRIP	= strip

CFLAGS = -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -pipe

ifeq "$(DEBUG)" ""
	CFLAGS += -O3
else
	CFLAGS += -O -g
endif

C?=0
ifeq "$(C)" "0"
check_source = 
else
CFLAGS+=-D_FORTIFY_SOURCE=2
check_source = echo "CPPCHECK $1..."; #cppcheck --quiet --std=c99 $1
endif

build/$(ARCH)/%.o: %.f95
	echo "CC $@..."
	mkdir -p $(dir $@)
	$(call check_source,$<)
	$(CF) $(CFLAGS) -c -o $@ $<

build/$(ARCH)/%.dep: %.f95
# Don't build dependencies if this is a clean build
ifneq ($(MAKECMDGOALS),clean)
	echo "DEP $<..."
	mkdir -p $(dir $@)
#	So, it turns out gfortran is bugged and that this doesn't generate
#	the dependencies like is expected. It works for gcc, but not for
#	gfortran. :(
	$(CF) $< $(CFLAGS) -cpp -MM -MG -MT $(basename $@).o -o $@
endif

# Turn off verbose output by default
V?=0
ifeq "$(V)" "0"
	MAKEFLAGS = --no-print-directory --quiet
endif
