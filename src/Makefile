# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = day2doy diff_apr diff_ret distance doy2day events extract hurricane island issifm jsec2time map_pert map_rad map_ret noise_pert noise_ret optimize_btd orbit overpass perturbation rayt ret2tab sampling select_ret spec2tab spec_ana spec_qual spec_synth sza time2jsec var1d variance volcano zm_ret

# Tests...
TESTS = pert_test rad_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include -I/usr/include -I/usr/include/hdf5/serial

# Library directories...
LIBDIR += -L ../libs/build/lib -L /usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial

# Static compilation...
STATIC ?= 0

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage testing...
COV ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)

# Compiler...
CC ?= gcc
MPICC ?= mpicc.openmpi

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums $(OPT) -g -fopenmp -mcmodel=medium

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm -lairs -lhdfeos -lmfhdf -ldf -lz -ljpeg -lnsl -lnetcdf -lm -ltirpc

# LDFLAGS for retrieval...
LDFLAGS_RET ?= $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lm

# Optimization information...
ifeq ($(INFO),1)
  CFLAGS += -fopt-info
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Compile for coverage testing...
ifeq ($(COV),1)
  CFLAGS += --coverage
endif

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen indent install lizard mkdocs strip uninstall

all: $(EXC) nlte retrieval
	rm -f *~

$(EXC): %: %.c libairs.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libairs.o jurassic.o $(LDFLAGS)

libairs.o: libairs.c libairs.h Makefile
	$(CC) $(CFLAGS) -c -o libairs.o libairs.c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

nlte: nlte.c jurassic.o
	$(MPICC) $(CFLAGS) -o nlte nlte.c jurassic.o $(LDFLAGS_RET)

retrieval: retrieval.c jurassic.o
	$(MPICC) $(CFLAGS) -o retrieval retrieval.c jurassic.o $(LDFLAGS_RET)

check: $(TESTS)

$(TESTS): all
	@(echo "\n===== Running \"$@\" ... =====") ; \
	  cd ../tests/$@ ; ./run.sh \
	  && (echo "\n===== Test \"$@\" passed! =====") \
	  || (echo "\n===== Test \"$@\" failed! =====" ; exit 1)

clean:
	rm -rf $(EXC) nlte retrieval *.o *.gcda *.gcno *.gcov coverage* *~

coverage:
	rm -f jurassic.gcda jurassic.gcno ; \
	lcov --capture --directory . --output-file=coverage.info ; \
	genhtml coverage.info --legend --frames --output-directory=coverage ; \
	gcov -abcfu *.c || true

cppcheck:
	cppcheck --enable=all ./

dist:
	cd ..  && tar cvfz airs_$(VERSION).tgz bin/* \
          docs/Doxyfile example/* libs/build.sh libs/*gz src/* tests/*

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC) nlte retrieval

uninstall:
	cd $(DESTDIR) && rm $(EXC)
