# Mirror-Neutron-Density-Evolution
This repo contains the Fortran implementation of Y. Kamyshkov, J. Ternullo, L. Varriano, and Z. Berezhiani's paper on neutron-mirror neutron oscillations in absorbing matter (https://doi.org/10.3390/sym14020230).

## Included External Libraries
1. uuid_module.f90 by @jacobwilliams (https://github.com/jacobwilliams/uuid-fortran)
2. f90getopt.F90 by @haniibrahim (https://github.com/haniibrahim/f90getopt)

All credit for these libraries go to their authors.

## Compiling the Code
This code is currently a WIP - it functions, but many things must still be done by hand (since it's not been worth my time to automate them). To compile, navigate to the desired program folder under `src` and check out the `Makefile`. The instructions are basically 
1. `make clean`
2. `make (call here)`

I've been compiling this with GCC's `gfortran-12.2.0`, but any version that supports Fortran 2003 *should* work without issue. 
The Cernlib library is *not* required.
