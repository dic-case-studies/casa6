FC = gfortran
OPT=-O2
FFLAGS = -fPIC $(OPT)

clean:
	rm -f *.o libblas.a libblas.so*

OBJS=	caxpy.o ccopy.o cdotc.o cdotu.o cgbmv.o cgemm.o cgemv.o \
	cgerc.o cgeru.o chbmv.o chemm.o chemv.o cher2.o cher2k.o \
	cher.o cherk.o chpmv.o chpr2.o chpr.o crotg.o cscal.o csrot.o \
	csscal.o cswap.o csymm.o csyr2k.o csyrk.o ctbmv.o ctbsv.o \
	ctpmv.o ctpsv.o ctrmm.o ctrmv.o ctrsm.o ctrsv.o dasum.o \
	daxpy.o dcabs1.o dcopy.o ddot.o dgbmv.o dgemm.o dgemv.o \
	dger.o dnrm2.o drot.o drotg.o drotm.o drotmg.o dsbmv.o \
	dscal.o dsdot.o dspmv.o dspr2.o dspr.o dswap.o dsymm.o \
	dsymv.o dsyr2.o dsyr2k.o dsyr.o dsyrk.o dtbmv.o dtbsv.o \
	dtpmv.o dtpsv.o dtrmm.o dtrmv.o dtrsm.o dtrsv.o dzasum.o \
	dznrm2.o icamax.o idamax.o isamax.o izamax.o lsame.o sasum.o \
	saxpy.o scabs1.o scasum.o scnrm2.o scopy.o sdot.o sdsdot.o \
	sgbmv.o sgemm.o sgemv.o sger.o snrm2.o srot.o srotg.o srotm.o \
	srotmg.o ssbmv.o sscal.o sspmv.o sspr2.o sspr.o sswap.o \
	ssymm.o ssymv.o ssyr2.o ssyr2k.o ssyr.o ssyrk.o stbmv.o \
	stbsv.o stpmv.o stpsv.o strmm.o strmv.o strsm.o strsv.o \
	xerbla_array.o xerbla.o zaxpy.o zcopy.o zdotc.o zdotu.o \
	zdrot.o zdscal.o zgbmv.o zgemm.o zgemv.o zgerc.o zgeru.o \
	zhbmv.o zhemm.o zhemv.o zher2.o zher2k.o zher.o zherk.o \
	zhpmv.o zhpr2.o zhpr.o zrotg.o zscal.o zswap.o zsymm.o \
	zsyr2k.o zsyrk.o ztbmv.o ztbsv.o ztpmv.o ztpsv.o ztrmm.o \
	ztrmv.o ztrsm.o ztrsv.o

static: $(OBJS)
	ar ruv libblas.a $(OBJS)
	ranlib libblas.a

shared: $(OBJS)
	cc -L$(PREFIX)/lib $(CFLAGS) -shared -Wl,-rpath,some-path-that-does-not-exist-for-the-purpose-of-reserving-a-buffer-for-future-adjustments-to-the-rpath -Wl,-soname,libblas-casatools.so.3 -o libblas-casatools.so.3.2.1 $(OBJS) -lm -lgfortran -lc
	cp libblas-casatools.so.3.2.1 $(PREFIX)/lib
	cd $(PREFIX)/lib && rm -f libblas-casatools.so.3 && ln -s libblas-casatools.so.3.2.1 libblas-casatools.so.3
	cd $(PREFIX)/lib && rm -f libblas.so && ln -s libblas-casatools.so.3.2.1 libblas.so
