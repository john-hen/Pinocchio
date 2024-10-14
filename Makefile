# Build configuration for GNU Make
#
# We will create three different builds in separate output folders.
#   release: Optimized for performance.
#   debug:   Includes symbols for interactive debugging.
#   memtest: Checks for memory leaks at run time.
#
# The latter is a separate build as it would otherwise lead to this error:
#   "LeakSanitizer does not work under ptrace (strace, gdb, etc)"
#
# Possibly add include directories and library options by calling like:
#   make INC="..." LIB="..."
#
# Or define $INC and $LIB as environment variables, possibly sourcing them
# from an "env" file for the site where the application is being deployed.


###########
# Options #
###########

# Activate, or not, compile-time options modifying program behavior.
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT
OPTIONS += -DPLC
OPTIONS += -DDO_NOT_REALLOC
#OPTIONS += -DROTATE_BOX
#OPTIONS += -DWHITENOISE
#OPTIONS += -DNO_RANDOM_MODULES
#OPTIONS += -DNORADIATION
#OPTIONS += -DSCALE_DEPENDENT_GROWTH

# We will create three different builds in separate output folders.
SOURCES         = $(notdir $(wildcard src/*.c))
OBJECTS         = $(patsubst %.c, %.o, $(SOURCES))
OBJECTS_RELEASE = $(addprefix build/release/objects/, $(OBJECTS))
OBJECTS_DEBUG   = $(addprefix build/debug/objects/,   $(OBJECTS))
OBJECTS_MEMTEST = $(addprefix build/memtest/objects/, $(OBJECTS))

# Use C compiler with MPI support by default.
CC = mpicc

# Define compiler flags for three different builds.
CFLAGS_RELEASE = -O3 -Wno-unused-result -Wno-format-overflow
CFLAGS_DEBUG   = -ggdb3 -Wall -fno-omit-frame-pointer
CFLAGS_MEMTEST = -fsanitize=address -g -Wall -fno-omit-frame-pointer

# Specify additional include directories, if any.
# The "override" keyword means we append to whatever value for INC
# the user may have passed in or defined in the shell environment.
override INC +=

# Specify the external libraries we depend on.
override LIB += -lfftw3_mpi -lfftw3 -lgslcblas -lgsl -lm

# Define install target using variable names as per the GNU coding standards.
# Installing Pinocchio is not necessary, as we can just call the executable
# from the "build" directory. Installation typically only makes sense when
# creating a lightweight container that contains nothing but the binaries.
DESTDIR =
PREFIX  = /usr
TARGET  = $(DESTDIR)$(PREFIX)


#########
# Rules #
#########

.PHONY: release debug memtest all install inspect clean

release: $(OBJECTS_RELEASE) build/release/pinocchio

debug:   $(OBJECTS_DEBUG)   build/debug/pinocchio

memtest: $(OBJECTS_MEMTEST) build/memtest/pinocchio

all: release debug memtest

install:
	install -d $(TARGET)/bin/
	install -m 777 build/release/pinocchio $(TARGET)/bin/

inspect:
	@echo "C compiler:      $(CC)"
	@echo "compile options: $(OPTIONS)"
	@echo "include options: $(INC)"
	@echo "library options: $(LIB)"
	@echo "install target:  $(TARGET)"

clean:
	rm -rf build/


build/release/objects/%.o: src/%.c src/*.h Makefile
	mkdir -p $(@D)
	$(CC) -c $(CFLAGS_RELEASE) $(INC) $(OPTIONS) -o $@ $<

build/release/pinocchio: $(OBJECTS_RELEASE) Makefile
	$(CC) $(CFLAGS_RELEASE) -o $@ $(OBJECTS_RELEASE) $(LIB)


build/debug/objects/%.o: src/%.c src/*.h Makefile
	mkdir -p $(@D)
	$(CC) -c $(CFLAGS_DEBUG) $(INC) $(OPTIONS) -o $@ $<

build/debug/pinocchio: $(OBJECTS_DEBUG) Makefile
	$(CC) $(CFLAGS_DEBUG) -o $@ $(OBJECTS_DEBUG) $(LIB)


build/memtest/objects/%.o: src/%.c src/*.h Makefile
	mkdir -p $(@D)
	$(CC) -c $(CFLAGS_MEMTEST) $(INC) $(OPTIONS) -o $@ $<

build/memtest/pinocchio: $(OBJECTS_MEMTEST) Makefile
	$(CC) $(CFLAGS_MEMTEST) -o $@ $(OBJECTS_MEMTEST) $(LIB)
