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

outA.txt logA.txt: mainA
	./$< 1> $@ 2> logA.txt

mainA: mainA.o qrMethods.o equal.o

mainA.o qrMethods.o: qrMethods.h


.PHONEY: B
B: outB.txt
	cat $<

outB.txt logB.txt: mainB
	./$< 1> $@ 2> logB.txt

mainB: mainB.o qrMethods.o equal.o

mainB.o qrMethods.o: qrMethods.h


.PHONEY: C
C: outC.txt
	cat $<

outC.txt logC.txt: mainC
	./$< 1> $@ 2> logC.txt

mainC: mainC.o qrMethods.o equal.o

mainC.o qrMethods.o: qrMethods.h


.PHONEY: clean
clean:
	$(RM) *.o *.txt mainA mainB mainC

