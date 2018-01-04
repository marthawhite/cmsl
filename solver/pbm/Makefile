# Make file for Proximal Bundle Method
# Copies the resulting mex into a secondary one,
# in case pbm calls itself.
# NOTE: mex files cannot call itself, or matlab crashes
# for nested optimization, pbm could call pbm, so instead
# in the code it calls pbm1

# Linux settings.
ifneq (,$(findstring Linux,$(OSTYPE)))
	FF = gfortran -O2 -fPIC
	MEX = /opt/matlab.R2008a/bin/mex
	FORTRAN_LIB = -lgfortran
	COPY_FILE = cp pbm.o pbm1.o
endif

# Mac OS X settings.
#ifneq (,$(findstring darwin,$(OSTYPE)))
 	FF = gfortran -O2 -fPIC
	MEX = mex
	FORTRAN_LIB = -L/usr/local/gfortran/lib -lgfortran
	COPY_FILE = cp pbm.mexmaci64 pbm1.mexmaci64
#endif

all:	pbm


pbmmex.o: pbmmex.f
	$(FF) -c pbmmex.f

mpbngc.o: mpbngc.f
	$(FF) -c mpbngc.f

plqdf1.o: plqdf1.f
	$(FF) -c plqdf1.f

pllpb2.o: pllpb2.f
	$(FF) -c pllpb2.f

pbm: pbmmex.o mpbngc.o plqdf1.o pllpb2.o
	$(MEX) -o pbm ${FORTRAN_LIB} pbm_driver.cpp objfunc.cpp pbmmex.o mpbngc.o plqdf1.o pllpb2.o

clean:
	rm -f pbm_driver.mex* *.o *.mexmaci64
