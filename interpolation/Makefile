CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: A B C


.PHONEY: A
A: plotLinterpRandom.svg plotLinterpSin.svg

plotLinterpRandom.svg plotLinterpSin.svg: plotLinterp.gpi dataA.dat
	gnuplot $<

dataA.dat: mainA
	./$< > $@

mainA: linterp.o


.PHONEY: B
B: plotQuad.svg

plotQuad.svg: plotQuad.gpi dataB.dat
	gnuplot $<

dataB.dat: mainB
	./$< > $@

mainB: mainB.o quadSpline.o

mainB.o quadSpline.o: quadSpline.h


.PHONEY: C
C: plotCubic.svg

plotCubic.svg: plotCubic.gpi dataC.dat
	gnuplot $<

dataC.dat: mainC
	./$< > $@

mainC: mainC.o cubicSpline.o

mainC.o cubicSpline.o: cubicSpline.h


.PHONEY: clean
clean:
	$(RM) *.o *.svg *.dat mainA mainB mainC
