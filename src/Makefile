# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = day2doy diff_apr diff_ret distance doy2day events extract hurricane island issifm map_pert map_rad map_ret noise_pert noise_ret optimize_btd orbit overpass perturbation rayt ret2tab sampling spec2tab spec_ana sza var1d variance volcano zm_ret

# Library directories...
LIBDIR = -L ../lib/build/lib

# Include directories...
INCDIR = -I ../lib/build/include

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc
MPICC = mpicc

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED $(DEFS) -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -g -mcmodel=medium

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm -lairs -lhdfeos -lGctp -lmfhdf -ldf -lz -ljpeg -lnsl -lnetcdf -lm

# LDFLAGS for retrieval...
LDFLAGS_RET = $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lm

# Profiling...
ifdef PROF
  CFLAGS += -O2 -pg 
else
  CFLAGS += -O3
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

dist:
	zip airs_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC) nlte retrieval
