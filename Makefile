API_VERSION = 2

OS = $(shell uname -s)

# I parse the version from the changelog. This version is generally something
# like 0.04-1 Here 0.04 is the main version and 1 is the debian package
# version. I only use the main version and strip leading 0s, so the above
# becomes 0.4
VERSION := $(shell sed -n 's/.*(\([0-9\.]*[0-9]\).*).*/\1/; s/\.0*/./g; p; q;' Changes)

ifeq ($(strip $(VERSION)),)
$(error "Couldn't parse version from debian/changelog")
endif

ifeq ($(OS),Darwin)
	CC = gcc # otherwise, CC defaults to clang
	SO_VERSION=$(API_VERSION).$(VERSION)
	TARGET_SO_BARE   = libdogleg.dylib
	TARGET_SO_FULL   = libdogleg.$(SO_VERSION).dylib
	TARGET_SO_SONAME = libdogleg.$(API_VERSION).dylib

# osx doesn't have DT_NEEDED or something, so I specify these explicitly
	LDLIBS += -lsuitesparseconfig -lamd -lcolamd -llapack -lblas
else
	SO_VERSION=$(API_VERSION).$(VERSION)
	TARGET_SO_BARE   = libdogleg.so
	TARGET_SO_FULL   = $(TARGET_SO_BARE).$(SO_VERSION)
	TARGET_SO_SONAME = $(TARGET_SO_BARE).$(API_VERSION)
	LDFLAGS += -Wl,-soname -Wl,libdogleg.so.$(API_VERSION)
endif

OPTFLAGS = -O3 -ffast-math -mtune=core2
FLAGS += -ggdb  -Wall -Wextra -MMD $(OPTFLAGS) -I/usr/include/suitesparse
CFLAGS += $(FLAGS) --std=gnu99

LDLIBS += -lcholmod
LDLIBS += -lm

HEADERS = dogleg.h

LIB_TARGETS = libdogleg.a $(TARGET_SO_BARE) $(TARGET_SO_FULL) $(TARGET_SO_SONAME)

MAN_SECTION = 3
MAN_TARGET = libdogleg.$(MAN_SECTION)

ALL_TARGETS = $(LIB_TARGETS) $(MAN_TARGET)

all: $(ALL_TARGETS)

%.a: dogleg.o
	ar rcvu $@ $^

$(TARGET_SO_FULL): dogleg-pic.o
	$(CC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

$(TARGET_SO_SONAME) $(TARGET_SO_BARE): $(TARGET_SO_FULL)
	ln -fs $^ $@

%-pic.o: %.c
	$(CC) -fPIC $(CFLAGS) -c -o $@ $<

$(MAN_TARGET): README.pod
	pod2man --center="libdogleg: Powell's dogleg method" --name=LIBDOGLEG --release="libdogleg $(VERSION)" --section=$(MAN_SECTION) $^ $@

ifdef DESTDIR
install: $(ALL_TARGETS)
	mkdir -p $(DESTDIR)/usr/lib/
	cp -P $(LIB_TARGETS) $(DESTDIR)/usr/lib/
	mkdir -p $(DESTDIR)/usr/include/
	install -m 0644 $(HEADERS) $(DESTDIR)/usr/include/
else
install:
	@echo "make install is here ONLY for the debian package. Do NOT run it yourself" && false
endif


sample: sample.o libdogleg.a
	$(CC) $^ $(LDLIBS) -o $@

sample.o: CFLAGS += -I.

clean:
	rm -f libdogleg.so* *.o *.a *.d $(ALL_TARGETS)

.PHONY: all clean install

-include *.d
