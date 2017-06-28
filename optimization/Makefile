CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B C


.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< > $@

mainA: mainA.o minimizationMethods.o qrMethods.o


.PHONEY: B
B: outB.txt plotB.svg
	cat $<

outB.txt: mainB mainBRoot
	./$< > $@
	@ echo "" >> $@
	./mainBRoot >> $@

mainB: mainB.o minimizationMethods.o qrMethods.o

mainBRoot: mainBRoot.o qrMethods.o newtonMethods.o newtonMethods.h

newtonMethods.o: newtonMethods.h qrMethods.h

plotB.svg: plotB.gpi outBFit.txt
	gnuplot $<

outBFit.txt: mainBFit
	./$< > $@

mainBFit: mainBFit.o minimizationMethods.o qrMethods.o


.PHONEY: C
C: outC.txt
	cat $<

outC.txt: mainC
	./$< > $@

mainC: mainC.o minimizationMethods.o qrMethods.o



minimizationMethods.o: qrMethods.o qrMethods.h

qrMethods.o: qrMethods.h



.PHONEY: clean
clean:
	$(RM) mainA mainB mainBRoot mainBFit mainC *.o *.txt *.svg
