CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: A B C


.PHONEY: A
A: plotLinterpRandom.svg plotLinterpSin.svg

plotLinterpRandom.svg plotLinterpSin.svg: plotLinterp.gpi dataA.dat
	gnuplot $<

dataA.dat: integLinterp
	./$< > $@

integLinterp: integLinterp.o linterp.o


.PHONEY: B
B: plotQuad.svg

plotQuad.svg: plotQuad.gpi dataB.dat
	gnuplot $<

dataB.dat: quadSplineTest
	./$< > $@

quadSplineTest: quadSplineTest.o quadSpline.o

quadSplineTest.o quadSpline.o: quadSpline.h


.PHONEY: C
C: plotCubic.svg

plotCubic.svg: plotCubic.gpi dataC.dat
	gnuplot $<

dataC.dat: cubicSplineTest
	./$< > $@

cubicSplineTest: cubicSplineTest.o cubicSpline.o

cubicSplineTest.o cubicSpline.o: cubicSpline.h


.PHONEY: clean
clean:
	$(RM) *.o *.svg *.dat integLinterp quadSplineTest cubicSplineTest
