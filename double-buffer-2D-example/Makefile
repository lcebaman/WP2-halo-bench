src = $(wildcard *.c)
obj = $(src:.c=.o)


CC = cc

BUILD :=overlap
$(info BUILD=${BUILD})

COMP :=yes
$(info COMP=${COMP})

LDFLAGS =

ifeq ($(BUILD),overlap)
	CFLAGS = -DOVERLAP
	EXE=overlap.exe
ifeq ($(COMP),no)
	CFLAGS += -DNOCOMPUTE
	EXE=overlap_nocomp.exe
endif

else
	CFLAGS=
	EXE=double.exe
ifeq ($(COMP),no)
	CFLAGS += -DNOCOMPUTE
	EXE=double_nocomp.exe
endif

endif	




%.o: %.c Makefile
	$(CC) $(CFLAGS) $< -c -o $@

${EXE}: $(obj)
	$(CC)$(LDFLAGS) -o $@ $^


.PHONY: clean

clean:
	rm -f $(obj) *.exe

