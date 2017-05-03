CFLAGS := -Wall -std=gnu99 $$(gsl-config --cflags) \
	  #-DHAVE_INLINE -Ofast -march=native \
	  #-DGSL_RANGE_CHECK_OFF \
	  #/proc/cpuinfo
LDLIBS := $$(gsl-config --libs)

.PHONEY: all
all: A B C


.PHONEY: A
A: outA.txt
	cat $<

outA.txt logA.txt: qrGsDecompTest
	./$< 1> $@ 2> logA.txt

qrGsDecompTest: qrGsDecompTest.o qrMethods.o equal.o

qrGsDecompTest.o qrMethods.o: qrMethods.h


.PHONEY: B
B: outB.txt
	cat $<

outB.txt logB.txt: qrGsInverseTest
	./$< 1> $@ 2> logB.txt

qrGsInverseTest: qrGsInverseTest.o qrMethods.o equal.o

qrGsInverseTest.o qrMethods.o: qrMethods.h


.PHONEY: C
C: outC.txt
	cat $<

outC.txt logC.txt: qrGivensTest
	./$< 1> $@ 2> logC.txt

qrGivensTest: qrGivensTest.o qrMethods.o equal.o

qrGivensTest.o qrMethods.o: qrMethods.h


.PHONEY: clean
clean:
	$(RM) *.o *.txt qrGsDecompTest qrGsInverseTest

