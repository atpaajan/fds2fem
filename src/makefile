# Makefile

INST_DIR=$(HOME)/bin

# Compilation instructions

FC= ifort
FFLAGS= -free -O3 -traceback
LDFLAGS=
OBJ= global.o string.o error.o misc.o math.o cfg_reader.o fds_dump.o fds_reader.o cfast_dump.o iso_reader.o \
		 fds_stats.o fem_dump.o fem_stats.o abaqus_reader.o abaqus_output.o \
		 ansys_reader.o ansys_output.o matching.o \
		 mapping.o mapping_iso.o main.o 

fds2fem_db : FFLAGS= -free -m64 -check -warn all -O0 -auto -WB -traceback -g -fpe0 -fltconsistency -no-wrap-margin

fds2fem: $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o fds2fem $(OBJ)

fds2fem_db: $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o fds2fem_db $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

install:
	cp fds2fem $(INST_DIR)

clean:
	rm -f *.o *.mod
