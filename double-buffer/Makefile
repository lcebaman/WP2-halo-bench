src = $(wildcard *.c)
obj = $(src:.c=.o)

LDFLAGS = 

myprog: $(obj)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(obj) myprog
