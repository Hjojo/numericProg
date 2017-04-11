CFLAGS := -Wall -std=gnu99
LDLIBS := -lm

.PHONEY: all
all: plotLinterp.svg

plotLinterp.svg plotLinterpSin.svg: plotLinterp.gpi data
	gnuplot $<

data: integLinterp
	./$< > $@

integLinterp: integLinterp.o linterp.o

.PHONEY: clean
clean:
	$(RM) *.o integLinterp data *.svg
