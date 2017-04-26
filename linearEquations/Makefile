CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B


.PHONEY: A
A: outA.txt
	cat $<

outA.txt logA.txt: qrGsDecompTest
	./$< 1> $@ 2> logA.txt

qrGsDecompTest: qrGsDecompTest.o qrMethods.o equal.o


.PHONEY: B
B: outB.txt
	cat $<

outB.txt logB.txt: qrGsInverseTest
	./$< 1> $@ 2> logB.txt

qrGsInverseTest: qrGsInverseTest.o qrMethods.o equal.o


.PHONEY: clean
clean:
	$(RM) *.o *.txt qrGsDecompTest qrGsInverseTest

