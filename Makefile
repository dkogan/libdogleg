include choose_mrbuild.mk
include $(MRBUILD_MK)/Makefile.common.header

PROJECT_NAME := dogleg

VERSION      := $(VERSION_FROM_PROJECT)
ABI_VERSION  := 2
TAIL_VERSION := $(VERSION)

LIB_SOURCES += \
  dogleg.c

BIN_SOURCES += \
  sample.c \
  test-misc.c

# I do not distribute ANY binaries
DIST_BIN_EXCEPT := *

DIST_INCLUDE += \
  dogleg.h

LDLIBS += -lcholmod -llapack
LDLIBS += -lm

CFLAGS += -Wall -Wextra -I/usr/include/suitesparse

sample.o: CFLAGS += -I.

check: sample test-misc
	./check.sh
.PHONY: check

MAN_SECTION := 3
MAN_TARGET  := libdogleg.$(MAN_SECTION)
DIST_MAN    := $(MAN_TARGET)

$(MAN_TARGET): README.pod
	pod2man --center="libdogleg: Powell's dogleg method" --name=LIBDOGLEG --release="libdogleg $(VERSION)" --section=$(MAN_SECTION) $^ $@

EXTRA_CLEAN += $(MAN_TARGET)

include $(MRBUILD_MK)/Makefile.common.footer
