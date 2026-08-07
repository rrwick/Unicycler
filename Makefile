# This makefile will build the C++ components of Unicycler.

# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the *.so file, which is required to run the aligner)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# Determine the platform
ifeq ($(shell uname), Darwin)
  PLATFORM = Mac
else
  PLATFORM = Linux
endif
$(info Platform: $(PLATFORM))


# Determine the compiler and version
COMPILER_HELP := $(shell $(CXX) --help | head -n 1)
ifneq (,$(findstring clang,$(COMPILER_HELP)))
    COMPILER = clang
else ifneq (,$(findstring g++,$(COMPILER_HELP)))
    COMPILER = g++
else ifneq (,$(findstring c++,$(COMPILER_HELP)))
    COMPILER = c++
else ifneq (,$(findstring Intel,$(COMPILER_HELP)))
    COMPILER = icpc
else
    COMPILER = unknown
endif
ifeq ($(COMPILER),clang)
  COMPILER_VERSION := $(shell $(CXX) --version | grep version | grep -o -m 1 "[0-9]\+\.[0-9]\+\.*[0-9]*" | head -n 1)
else
  COMPILER_VERSION := $(shell $(CXX) -dumpfullversion 2> /dev/null)
  ifeq ($(COMPILER_VERSION),)
      $(info Falling back to -dumpversion as compiler did not support -dumpfullversion)
      COMPILER_VERSION := $(shell $(CXX) -dumpversion)
  endif
endif
COMPILER_VERSION_NUMBER := $(shell echo $(COMPILER_VERSION) | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/')
$(info Compiler: $(COMPILER) $(COMPILER_VERSION))


# Clang versions earlier than 3.5 won't work with Seqan.
ifeq ($(COMPILER),clang)
  CLANG_340_OR_MORE := $(shell expr $(COMPILER_VERSION_NUMBER) \>= 30500)
  ifeq ($(CLANG_340_OR_MORE),0)
    $(error Unicycler requires Clang 3.5 or greater)
  endif
endif


# GCC versions earlier than 4.9.1 won't work with Seqan.
ifeq ($(COMPILER),g++)
  GCC_491_OR_MORE := $(shell expr $(COMPILER_VERSION_NUMBER) \>= 40901)
  ifeq ($(GCC_491_OR_MORE),0)
    $(error Unicycler requires GCC version 4.9.1 or greater)
  endif
endif


# CXXFLAGS can be overridden by the user.
CXXFLAGS    ?= -Wall -Wextra -pedantic -mtune=native


# These flags are required for the build to work.
FLAGS        = -std=c++14 -Iunicycler/include -fPIC
LDFLAGS      += -shared -lz


# Platform-specific stuff (for Seqan)
ifeq ($(PLATFORM), Linux)
  FLAGS     += -lrt -lpthread
endif


# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -DNDEBUG


TARGET       = unicycler/cpp_functions.so
SHELL        = /bin/sh
SOURCES      = $(shell find unicycler -name "*.cpp")
HEADERS      = $(shell find unicycler -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)


# Linux needs '-soname' while Mac needs '-install_name'
ifeq ($(PLATFORM), Mac)
  SONAME       = -install_name
else  # Linux
  SONAME       = -soname
endif


.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)


.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
