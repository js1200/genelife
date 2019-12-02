# makefile for making C library genelife_lib.so
#
# written by John S. McCaskill Nov 2019#
#
# use make -f Makefile_xterm to make standalone C xterm application with primitive graphics without python
#
TARGET_LIB ?= libgenelife

BUILD_DIR ?= ./buildlib
SRC_DIRS ?= ./genelifec ./genelifeclib

OS_LDFLAG :=
OS_NAME := $(shell uname -s)
ifeq ($(OS_NAME),Linux)
	LDLIBS := $(TARGET_LIB).so
endif
ifeq ($(OS_NAME),Darwin)
	LDLIBS := $(TARGET_LIB).dylib
	OS_LDFLAG = -dynamiclib
endif

SRCS := $(shell find $(SRC_DIRS) -name *.c)

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CFLAGS ?= $(INC_FLAGS) -fPIC -g -Wall -std=gnu99 -MMD -MP # C flags,    add -Wextra to see mixed integer type warnings etc
LDFLAGS := $(OS_LDFLAG) -lm -shared   # linking flags

$(LDLIBS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) $(LDLIBS)

-include $(DEPS)

MKDIR_P ?= mkdir -p
