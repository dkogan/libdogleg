API_VERSION = 1
VERSION := $(shell dpkg-parsechangelog | awk '/^Version/{ gsub("-.*","",$$2); print $$2}' | perl -pe 's/\.0*([1-9])/.$$1/g')

ifeq ($(strip $(VERSION)),)
$(error "Couldn't parse version from debian/changelog")
endif

SO_VERSION=$(API_VERSION).$(VERSION)

OPTFLAGS = -O3 -ffast-math -mtune=core2
FLAGS += -ggdb  -Wall -Wextra -MMD $(OPTFLAGS) -I/usr/include/suitesparse
CFLAGS = $(FLAGS) --std=gnu99

LDFLAGS  += -ggdb

LDLIBS += -lcholmod
LDLIBS += -lm

HEADERS = dogleg.h
TARGET_SO = libdogleg.so.$(SO_VERSION)
LIB_TARGETS = libdogleg.a $(TARGET_SO)

MAN_SECTION = 3
MAN_TARGET = libdogleg.$(MAN_SECTION)

ALL_TARGETS = $(LIB_TARGETS) $(MAN_TARGET)

all: $(ALL_TARGETS)

%.a: dogleg.o
	ar rcvu $@ $^

%.so.$(SO_VERSION): dogleg-pic.o
	$(CC) $(LDLIBS) -shared  $^ -Wl,-soname -Wl,libdogleg.so.$(API_VERSION) -o $@

%-pic.o: %.c
	$(CC) -fPIC $(CFLAGS) -c -o $@ $<

$(MAN_TARGET): README.pod
	pod2man --center="libdogleg: Powell's dogleg method" --name=LIBDOGLEG --release="libdogleg $(VERSION)" --section=$(MAN_SECTION) $^ $@

ifdef DESTDIR
install:
	mkdir -p $(DESTDIR)/usr/lib/
	install -m 0644 $(LIB_TARGETS) $(DESTDIR)/usr/lib/
	cd $(DESTDIR)/usr/lib/ && \
	ln -fs $(TARGET_SO) libdogleg.so.$(API_VERSION) && \
	ln -fs $(TARGET_SO) libdogleg.so && \
	cd -
	mkdir -p $(DESTDIR)/usr/include/
	install -m 0644 $(HEADERS) $(DESTDIR)/usr/include/
else
install:
	@echo "make install is here ONLY for the debian package. Do NOT run it yourself" && false
endif


clean:
	rm -f libdogleg.so* *.o *.a *.d $(ALL_TARGETS)

.PHONY: all clean install

-include *.d
