CC = gfortran
CFLAGS = -O3 -cpp -Wall -Wextra -pedantic -fcheck=all -fbacktrace -Wno-unused-variable -Wno-unused-function -fdefault-real-8 -fno-align-commons
#CERNLPATH=/usr/lib/cernlib/2006-g77/lib
CERNLPATH=/usr/lib64/

mathlib = -lm
cernlib = $(CERNLPATH)/libpacklib.a $(CERNLPATH)/libkernlib.a $(CERNLPATH)/libmathlib.a
#cernlib = libpacklib.a libkernlib.a libmathlib.a

libpath= ../lib
libobj = $(libpath)/$(wildcard $(libpath)/*.o)
libmod = $(libpath)/$(wildcard $(libpath)/*.mod)

banforb4c: banforb4c.o
	#$(CC) banforb4c.o $(cernlib) -I$(libpath) $(libobj) -o ~/banforb4c -static
	$(CC) banforb4c.o -I$(libpath) $(libobj) -o ~/banforb4c -static

banforb4c.o: banforb4c.f
	#$(CC) $(CFLAGS) $(CERNFLAGS) -I$(libpath) -c banforb4c.f -static
	$(CC) $(CFLAGS) -I$(libpath) -c banforb4c.f -static

banformap: banformap.o
	$(CC) banformap.o $(cernlib) -o banformap -static

banformap.o: banformap.f
	$(CC) $(CFLAGS) $(CERNFLAGS) -c banformap.f -static

banfor: banfor.o
	$(CC) banfor.o $(cernlib) -o banfor -static

banfor.o: banfor.f
	$(CC) $(CFLAGS) $(CERNFLAGS) -c banfor.f -static

clean:
	rm -f *.o core.* log.txt
