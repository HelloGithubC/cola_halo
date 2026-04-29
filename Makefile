#
# cola_halo
#   

# Define OPENMP to enable MPI+OpenMP hybrid parallelization
OPENMP  = -fopenmp # -fopenmp for gcc
# Note that MacOS llvm does not support OPENMP

#CC      = mpicc -std=c99 
CC      = mpicc -std=c99 ## mpicc -std=c99 
WOPT    ?= -Wall
CFLAGS  := -O3 $(WOPT) $(OPENMP) -Wall
LIBS    := -lm

# Source and object directories
SRCDIR  = src
OUTDIR  = output

# Compile options
CFLAGS += -DUSE_LUA   # use LUA programming language for the parameter file

# Define paths of FFTW3 & GSL libraries if necessary.

LUA_DIR   ?= /usr/include/lua5.4 #e.g. /opt/local
FFTW3_DIR ?= #/public1/home/scb4289/software/fftw-3.3.9 #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #/public1/soft/gsl/2.5 #e.g. /Users/jkoda/Research/opt/gcc/gsl

CFLAGS += -I$(LUA_DIR) -I$(SRCDIR)

EXEC = cola_halo # halo
all: $(EXEC)

# Source files
SRCS := $(wildcard $(SRCDIR)/*.c)

# Object files with paths
OBJS := $(OUTDIR)/main.o
OBJS += $(OUTDIR)/lpt.o $(OUTDIR)/msg.o $(OUTDIR)/power.o $(OUTDIR)/confirm_param.o
OBJS += $(OUTDIR)/solve_growth.o
OBJS += $(OUTDIR)/pm.o $(OUTDIR)/cola.o $(OUTDIR)/fof.o $(OUTDIR)/comm.o $(OUTDIR)/move.o $(OUTDIR)/move_min.o
OBJS += $(OUTDIR)/write.o $(OUTDIR)/timer.o $(OUTDIR)/mem.o
OBJS += $(OUTDIR)/subsample.o $(OUTDIR)/coarse_grid.o

LIBS += -lgsl -lgslcblas
LIBS += -lfftw3f_mpi -lfftw3f ${OPENMP}
#LIBS += -L/public1/home/scb4289/software/fftw3_float_openmpi/lib/fftw3f_mpi.a -L/public1/home/scb4289/software/fftw3_float_openmpi/lib/fftw3f.a -L/public1/home/scb4289/software/fftw3_float_openmpi/lib

ifeq (,$(findstring -DUSE_LUA, $(CFLAGS)))
  OBJS += $(OUTDIR)/read_param.o
else
  OBJS += $(OUTDIR)/read_param_lua.o
  LIBS += -llua -ldl
endif

ifdef OPENMP
  LIBS += -lfftw3f_omp
  #LIBS += -L/public1/home/scb4289/software/fftw-3.3.9/lib/fftw3f_omp.a
  #LIBS += -L/public1/home/scb4289/software/fftw-3.3.9/lib
  #LIBS += -lfftw3f_threads       # for thread parallelization instead of omp
endif

cola_halo: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $@

# Pattern rule for compiling .c to .o
$(OUTDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Special rule for move_min.c and move_min.h generation
$(SRCDIR)/move_min.c: $(SRCDIR)/move.c
	echo "// This code is automatically generated from $<" > $@
	cat $< | sed -e 's/Particles/Snapshot/g' -e 's/Particle/ParticleMinimum/g' -e 's/move_particles2/move_particles2_min/' >> $@

$(SRCDIR)/move_min.h: $(SRCDIR)/move.h
	echo "// This code is automatically generated from $<" > $@
	cat $< | sed -e 's/Particles/Snapshot/g' -e 's/Particle/ParticleMinimum/g' -e 's/move_particles2/move_particles2_min/' >> $@

# Dependency rules (updated paths)
$(OUTDIR)/main.o: $(SRCDIR)/main.c $(SRCDIR)/parameters.h $(SRCDIR)/lpt.h $(SRCDIR)/particle.h $(SRCDIR)/msg.h $(SRCDIR)/power.h $(SRCDIR)/comm.h $(SRCDIR)/pm.h \
  $(SRCDIR)/cola.h $(SRCDIR)/fof.h $(SRCDIR)/write.h $(SRCDIR)/timer.h $(SRCDIR)/mem.h $(SRCDIR)/move.h $(SRCDIR)/subsample.h $(SRCDIR)/coarse_grid.h $(SRCDIR)/solve_growth.h
$(OUTDIR)/solve_growth.o: $(SRCDIR)/solve_growth.c $(SRCDIR)/particle.h $(SRCDIR)/parameters.h $(SRCDIR)/solve_growth.h
$(OUTDIR)/cola.o: $(SRCDIR)/cola.c $(SRCDIR)/particle.h $(SRCDIR)/msg.h $(SRCDIR)/cola.h $(SRCDIR)/timer.h $(SRCDIR)/solve_growth.h
$(OUTDIR)/comm.o: $(SRCDIR)/comm.c $(SRCDIR)/msg.h $(SRCDIR)/comm.h
$(OUTDIR)/fof.o: $(SRCDIR)/fof.c $(SRCDIR)/particle.h $(SRCDIR)/msg.h $(SRCDIR)/comm.h $(SRCDIR)/timer.h $(SRCDIR)/move_min.h
$(OUTDIR)/lpt.o: $(SRCDIR)/lpt.c $(SRCDIR)/msg.h $(SRCDIR)/power.h $(SRCDIR)/particle.h
$(OUTDIR)/mem.o: $(SRCDIR)/mem.c $(SRCDIR)/fof.h $(SRCDIR)/particle.h $(SRCDIR)/mem.h $(SRCDIR)/msg.h $(SRCDIR)/comm.h
$(OUTDIR)/move.o: $(SRCDIR)/move.c $(SRCDIR)/msg.h $(SRCDIR)/move.h $(SRCDIR)/particle.h $(SRCDIR)/comm.h
$(OUTDIR)/msg.o: $(SRCDIR)/msg.c $(SRCDIR)/msg.h
$(OUTDIR)/pm.o: $(SRCDIR)/pm.c $(SRCDIR)/pm.h $(SRCDIR)/particle.h $(SRCDIR)/msg.h $(SRCDIR)/comm.h $(SRCDIR)/timer.h
$(OUTDIR)/power.o: $(SRCDIR)/power.c $(SRCDIR)/msg.h
$(OUTDIR)/read_param_lua.o: $(SRCDIR)/read_param_lua.c $(SRCDIR)/parameters.h $(SRCDIR)/msg.h $(SRCDIR)/solve_growth.h 
$(OUTDIR)/timer.o: $(SRCDIR)/timer.c $(SRCDIR)/msg.h $(SRCDIR)/timer.h
$(OUTDIR)/write.o: $(SRCDIR)/write.c $(SRCDIR)/msg.h $(SRCDIR)/comm.h $(SRCDIR)/write.h $(SRCDIR)/particle.h

$(OUTDIR)/move_min.o: $(SRCDIR)/move_min.c $(SRCDIR)/msg.h $(SRCDIR)/move_min.h $(SRCDIR)/particle.h $(SRCDIR)/comm.h

#
# "halo" -- cola_halo without cola, only does FoF etc.
#
OBJS2 := $(OUTDIR)/halo_main.o $(OUTDIR)/read.o
OBJS2 += $(OUTDIR)/read_param_lua.o $(OUTDIR)/msg.o $(OUTDIR)/fof.o $(OUTDIR)/comm.o $(OUTDIR)/move_min.o
OBJS2 += $(OUTDIR)/write.o $(OUTDIR)/timer.o $(OUTDIR)/mem.o 
OBJS2 += $(OUTDIR)/subsample.o $(OUTDIR)/coarse_grid.o

halo: $(OBJS2)
	$(CC) $(OBJS2) $(LIBS) -o $@

.PHONY: clean run dependence

clean :
	/bin/rm -f $(EXEC) $(OBJS) $(OBJS2) $(SRCDIR)/move_min.c $(SRCDIR)/move_min.h

run:
	mpirun -n 2 ./cola_halo param.lua

dependence:
	gcc -MM -MG $(SRCDIR)/*.c
