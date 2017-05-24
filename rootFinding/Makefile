CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B


.PHONEY: A
A: outA.txt
	cat $<

outA.txt: mainA
	./$< > $@

mainA: mainA.o newtonMethods.o qrMethods.o

mainA.o: newtonMethods.h

newtonMethods.o: newtonMethods.h qrMethods.h

qrMethods.o: qrMethods.h


.PHONEY: B
B: outB.txt
	cat $<

outB.txt: mainB
	./$< > $@

mainB: mainB.o newtonMethods.o qrMethods.o


.PHONEY: clean
clean:
	$(RM) mainA mainB *.o *.txt *.svg
