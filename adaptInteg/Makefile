CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: A B


.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< > $@

mainA: mainA.o adaptInteg.o

mainA.o: adaptInteg.h


.PHONEY: B
B: outB.txt
	cat $<

outB.txt: mainB
	./$< > $@

mainB: mainB.o adaptInteg.o

mainB.o: adaptInteg.h


adaptInteg.o: adaptInteg.h


.PHONEY: clean
clean:
	$(RM) *.o *.txt mainA mainB
