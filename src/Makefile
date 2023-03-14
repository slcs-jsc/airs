# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = day2doy diff_apr diff_ret distance doy2day events extract hurricane island issifm map_pert map_rad map_ret noise_pert noise_ret optimize_btd orbit overpass perturbation rayt ret2tab sampling select_ret spec2tab spec_ana spec_synth sza var1d variance volcano zm_ret

# Installation directory...
DESTDIR ?= ../bin

# Library directories...
LIBDIR += -L ../lib/build/lib

# Include directories...
INCDIR += -I ../lib/build/include

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Compiler...
CC ?= gcc
MPICC ?= mpicc

# CFLAGS...
CFLAGS ?= $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums $(OPT) -g -fopenmp -mcmodel=medium

# LDFLAGS...
LDFLAGS ?= $(LIBDIR) -lgsl -lgslcblas -lm -lairs -lhdfeos -lGctp -lmfhdf -ldf -lz -ljpeg -lnsl -lnetcdf -lm -ltirpc

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

# Static compilation...
ifeq ($(STATIC),1)
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC) nlte retrieval
	rm -f *~

$(EXC): %: %.c libairs.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libairs.o jurassic.o $(LDFLAGS) -static

libairs.o: libairs.c libairs.h Makefile
	$(CC) $(CFLAGS) -c -o libairs.o libairs.c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

nlte: nlte.c jurassic.o
	$(MPICC) $(CFLAGS) -o nlte nlte.c jurassic.o $(LDFLAGS_RET)

retrieval: retrieval.c jurassic.o
	$(MPICC) $(CFLAGS) -o retrieval retrieval.c jurassic.o $(LDFLAGS_RET)

clean:
	rm -f $(EXC) nlte retrieval *.o *~

check_cppcheck:
	cppcheck --enable=all ./

check_lizard:
	lizard -s cyclomatic_complexity

dist:
	zip airs_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC) nlte retrieval
