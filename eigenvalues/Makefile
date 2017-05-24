CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B C


.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< 12 > $@

mainA: mainA.o eigenvalueMethods.o equal.o


.PHONEY: B
B: outB.txt
	cat $<

outB.txt: mainB mainBCompare
	./$< 12 4 > $@
	./mainBCompare 50 >> $@

mainB: mainB.o eigenvalueMethods.o

mainBCompare: mainBCompare.o eigenvalueMethods.o


.PHONEY: C
C: outC.txt
	cat $<

outC.txt: mainC Makefile
	./$< 100 1 > $@
	\time -f "Took %U seconds" -ao $@ ./$< 100 1 > /dev/null
	./$< 100 2 >> $@
	\time -f "Took %U seconds" -ao $@ ./$< 100 2 > /dev/null

mainC: mainC.o eigenvalueMethods.o


.PHONEY: clean
clean:
	$(RM) mainA *.o *.txt mainB mainBCompare mainC
