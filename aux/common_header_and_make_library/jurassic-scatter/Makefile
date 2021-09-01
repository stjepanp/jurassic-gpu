
#gcc prog.c vector.c mul_vector.c -L /p/software/jusuf/stages/2020/modules/all/Compiler/GCC/9.3.0 -lstdc++ -L /p/software/jusuf/stages/2020/software/CUDA/11.0/lib64 -lcudart  -L ../jurassic-gpu -lutil -I../jurassic-gpu/inter_folder -I ../common_header -o prog
# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = prog

# Library directories...
#LIBDIR = -L /usr/lib64 
LIBDIR = 
LIBDIR += -L /p/software/jusuf/stages/2020/modules/all/Compiler/GCC/9.3.0
LIBDIR += -L /p/software/jusuf/stages/2020/software/CUDA/11.0/lib64
LIBDIR += -L ../jurassic-gpu

# Include directories...
INCDIR = 
#INCDIR = -I /usr/local/include
INCDIR += -I ../common_header
INCDIR += -I ../jurassic-gpu/inter_folder

# OpenMP...
OMP = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler with/without openMP and without MPI...
# -----------------------------------------------------------------------------

# Compiler...
CC = gcc

# CFLAGS...
CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE

# Flags for OpenMP...
ifdef OMP
CFLAGS += -fopenmp
endif

# -----------------------------------------------------------------------------
# Set flags for open MPI C wrapper compiler with/without openMP...
# -----------------------------------------------------------------------------
# # Compiler...
# CC = mpicc 

# # CFLAGS...
# CFLAGS = $(INCDIR) -pedantic-errors -Wall -Wextra -Werror -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Wno-long-long -Winline -fno-common -DHAVE_INLINE

# # Flags for OpenMP...
# ifdef OMP
# CFLAGS += -fopenmp
# endif
# -----------------------------------------------------------------------------

CFLAGS += -D WORK_QUEUE


# Flags for profiling...
ifdef PROF
CFLAGS += -O2 -g -pg
else
CFLAGS += -O3
endif

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm -lstdc++ -lcudart -lutil

OBJECTS = vector.c mul_vector.c

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c $(OBJECTS) ../jurassic-gpu/libutil.a
	$(CC) $(CFLAGS) -o $@ $@.c $(OBJECTS) $(LDFLAGS)

../jurassic-gpu/libutil.a:
	make libutil.a -C ../jurassic-gpu

# Generic rules for all object files...
%.o:%.c %.d Makefile
	$(CC) $(CFLAGS) -c $<

bak:
	tar czf ../bak/jurassic_`date +"%y%m%d%H%M"`.tgz Makefile *.c *.h

clean:
	rm -f $(EXC) *.o *~
