CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B C


.PHONEY: A
A: plotA.svg

plotA.svg: plotA.gpi outA.dat
	gnuplot $<

outA.dat: mainA
	./$< > $@

mainA: mainA.o leastSquareMethods.o qrMethods.o

mainA.o: leastSquareMethods.h

leastSquareMethods.o: leastSquareMethods.h qrMethods.h

qrMethods.o: qrMethods.h


.PHONEY: B
B: plotB.svg

plotB.svg: plotB.gpi outB.dat
	gnuplot $<

outB.dat: mainB
	./$< > $@

mainB: mainB.o leastSquareErrorMethods.o jacobiEigenvalue.o qrMethods.o

mainB.o: leastSquareErrorMethods.h

leastSquareErrorMethods.o: leastSquareErrorMethods.h qrMethods.h jacobiEigenvalue.o

qrMethods.o: qrMethods.h


.PHONEY: C
C: plotC.svg

plotC.svg: plotC.gpi outC.dat
	gnuplot $<

outC.dat: mainC
	./$< > $@

mainC: mainC.o leastSquareErrorMethods.o jacobiEigenvalue.o qrMethods.o

mainC.o: leastSquareErrorMethods.h

leastSquareErrorMethods.o: leastSquareErrorMethods.h qrMethods.h jacobiEigenvalue.o

qrMethods.o: qrMethods.h


.PHONEY: clean
clean:
	$(RM) mainA mainB mainC *.o *.dat *.svg
