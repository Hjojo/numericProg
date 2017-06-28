CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B


.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< > $@

mainA: mainA.o ODEMethods.o

mainA.o: ODEMethods.h


.PHONEY: B
B: outB.txt
	cat $<

outB.txt: mainB
	./$< > $@

mainB: mainB.o ODEMethods.o

mainB.o: ODEMethods.h


ODEMethods.o: ODEMethods.h


.PHONEY: clean
clean:
	$(RM) *.txt *.o mainA mainB
