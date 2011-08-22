SO_VERSION=1.1.1

OPTFLAGS = -O3 -ffast-math -mtune=core2
FLAGS += -ggdb  -Wall -Wextra -MMD $(OPTFLAGS) -I/usr/include/suitesparse
CFLAGS = $(FLAGS) --std=gnu99

LDFLAGS  += -ggdb

LDLIBS += -lcholmod
LDLIBS += -lm

HEADERS = dogleg.h
TARGET_SO = libdogleg.so.$(SO_VERSION)
TARGETS = libdogleg.a $(TARGET_SO)

all: $(TARGETS)

%.a: dogleg.o
	ar rcvu $@ $^

%.so.$(SO_VERSION): dogleg-pic.o
	$(CC) $(LDLIBS) -shared  $^ -Wl,-soname -Wl,libdogleg.so.1 -o $@

%-pic.o: %.c
	$(CC) -fPIC $(CFLAGS) -c -o $@ $<

ifdef DESTDIR
install:
	mkdir -p $(DESTDIR)/usr/lib/
	install -m 0644 -s $(TARGETS) $(DESTDIR)/usr/lib/
	cd $(DESTDIR)/usr/lib/
	ln -fs $(TARGET_SO) libdogleg.so.1
	ln -fs $(TARGET_SO) libdogleg.so
	cd -

	mkdir -p $(DESTDIR)/usr/include/
	install -m 0644 $(HEADERS) $(DESTDIR)/usr/include/
else
install:
	@echo "make install is here ONLY for the debian package. Do NOT run it yourself" && false
endif


clean:
	rm -f libdogleg.so* *.o *.a *.d

.PHONY: all clean install

-include *.d
