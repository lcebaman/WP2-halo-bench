src = $(wildcard *.c)
obj = $(src:.c=.o)

CC = mpicc
LDFLAGS = 

myprog.exe: $(obj)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(obj) myprog.exe
