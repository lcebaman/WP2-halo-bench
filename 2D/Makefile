src = $(wildcard *.c)
obj = $(src:.c=.o)

CC = cc
LDFLAGS = 

myprog.exe: $(obj)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(obj) myprog.exe
