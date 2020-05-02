API_VERSION = 2

OS = $(shell uname -s)

# I parse the version from the changelog. This version is generally something
# like 0.04 .I strip leading 0s, so the above becomes 0.4
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

# Add my optimization flags if none already exist. These could exist if debian
# packages are built.
OPTFLAGS := $(if $(filter -O%, $(CFLAGS)),,-O3 -ffast-math -mtune=core2) -fno-omit-frame-pointer

FLAGS += -ggdb  -Wall -Wextra -MMD $(OPTFLAGS) -I/usr/include/suitesparse
CFLAGS += $(FLAGS) --std=gnu99

LDLIBS += -lcholmod -lsuitesparseconfig -llapack
LDLIBS += -lm

LIB_TARGETS = libdogleg.a $(TARGET_SO_BARE) $(TARGET_SO_FULL) $(TARGET_SO_SONAME)

MAN_SECTION = 3
MAN_TARGET = libdogleg.$(MAN_SECTION)

ALL_TARGETS = $(LIB_TARGETS) $(MAN_TARGET) sample

all: $(ALL_TARGETS)

libdogleg.a: dogleg.o
	ar rcvu $@ $^

$(TARGET_SO_FULL): dogleg-pic.o
	$(CC) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

$(TARGET_SO_SONAME) $(TARGET_SO_BARE): $(TARGET_SO_FULL)
	ln -fs $^ $@

%-pic.o: %.c
	$(CC) $(CPPFLAGS) -fPIC $(CFLAGS) -c -o $@ $<

$(MAN_TARGET): README.pod
	pod2man --center="libdogleg: Powell's dogleg method" --name=LIBDOGLEG --release="libdogleg $(VERSION)" --section=$(MAN_SECTION) $^ $@


sample: sample.o libdogleg.a
	$(CC) $^ $(LDLIBS) -o $@

sample.o: CFLAGS += -I.

clean:
	rm -f libdogleg.so* *.o *.a *.d $(ALL_TARGETS)

.PHONY: all clean

-include *.d
