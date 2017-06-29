CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: plotRND.svg plotSin.svg


plotRND.svg: plotRND.gpi outRND.txt
	gnuplot $<

outRND.txt: mainRND
	./$< > $@

mainRND: mainRND.o cubicSpline.o cubicSubSpline.o quadSpline.o

mainRND.o: cubicSpline.h cubicSubSpline.h quadSpline.h



plotSin.svg: plotSin.gpi outSin.txt
	gnuplot $<

outSin.txt: mainSin
	./$< > $@

mainSin: mainSin.o cubicSpline.o cubicSubSpline.o quadSpline.o

mainSin.o: cubicSpline.h cubicSubSpline.h quadSpline.h



cubicSpline.o: cubicSpline.h

cubicSubSpline.o: cubicSpline.h quadSpline.h

quadSpline.o: quadSpline.h


.PHONEY:
clean:
	$(RM) *.o *.txt *.svg mainRND mainSin
