CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: A B



.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< > $@

mainA: mainA.o plainMonteCarlo.o


.PHONEY: B
B: plotB.svg

plotB.svg: plotB.gpi outB.txt
	gnuplot $<

outB.txt: mainB
	./$< > $@

mainB: mainB.o plainMonteCarlo.o


.PHONEY: clean
clean:
	$(RM) *.txt *.o *.svg mainA mainB
