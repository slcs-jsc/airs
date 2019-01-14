# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = arpp autocorr bands cfc collect day2doy distance doy2day erafm events extract get_wave_pert get_wave_synth hurricane island map_pert map_rad map_ret noise_pert noise_ret optimize_si orbit overpass pca perturbation rayt ret2tab sampling so2 spec2tab statistics sza umfm var1d var3d variance volcano wrffm zm_ret

# Library directories...
LIBDIR = -L ../lib/build/lib -L $(HOME)/lib

# Include directories...
INCDIR = -I ../lib/build/include -I $(HOME)/include

# Linking...
STATIC = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc
MPICC = mpicc

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED $(DEFS) -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm -lairs -lhdfeos -lGctp -lmfhdf -ldf -lz -ljpeg -lnsl -lnetcdf -lm

# LDFLAGS for retrieval...
LDFLAGS_RET = $(LIBDIR) -lgsl -lgslcblas -lm -lnetcdf -lm

# Profiling...
ifdef PROF
  CFLAGS += -O2 -g -pg 
else
  CFLAGS += -O3
endif

# Linking...
ifdef STATIC
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC) retrieval
	rm -f *~

$(EXC): %: %.c libairs.o jurassic.o
	$(CC) $(CFLAGS) -o $@ $< libairs.o jurassic.o $(LDFLAGS)

libairs.o: libairs.c libairs.h Makefile
	$(CC) $(CFLAGS) -c -o libairs.o libairs.c

jurassic.o: jurassic.c jurassic.h Makefile
	$(CC) $(CFLAGS) -c -o jurassic.o jurassic.c

bak:
	tar chzf airs_`date +"%y%m%d%H%M"`.tgz Doxyfile Makefile *.c *.h

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

clean:
	rm -f $(EXC) *.o *~

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

retrieval: retrieval.c jurassic.o
	$(MPICC) $(CFLAGS) -o retrieval retrieval.c jurassic.o $(LDFLAGS_RET)
