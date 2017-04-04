CFLAGS:= -Wall -std=gnu99

.PHONEY: all
all: A B C

.PHONEY: A
A: out.A.txt
	cat $<

out.A.txt: helloUser
	./$< > $@

.PHONEY: B
B: out.B.txt
	cat $<

out.B.txt: helloUserFun
	./$< > $@

helloUserFun: helloUserFun.o printHello.o printMe.o

.PHONEY: C
C: check

.PHONEY: check
check: out.A.txt out.B.txt correct.txt
	@echo "Testing if target A is correct..."
	@diff --brief out.A.txt correct.txt
	@echo "Target A is fine!"
	@echo "Testing if target B is correct..."
	@diff --brief out.B.txt correct.txt
	@echo "Target B is fine!"

correct.txt:
	echo "Hello $$(whoami)" > correct.txt


.PHONEY: clean
clean:
	$(RM) *.txt helloUser *.o helloUserFun
