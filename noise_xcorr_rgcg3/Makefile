
CC=g++
FC=gfortran
FFLAGS=-g
CFLAGS=-Wall -g
LIBBOOST = -lboost_filesystem -lboost_program_options -lboost_system
INCBOOST = -I /usr/include/boost
LGSL = -lgsl -lgslcblas
LFFTW = -lfftw3

BINDIR = $(HOME)/software/bin
INSTALL = install

PROGS = split_data noise_preprocess noise_xcorr

.SUFFIXES:
.SUFFIXES: .o .cpp .f

.cpp.o :
	$(CC) $(CFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

all: $(PROGS)

split_data: split_data.cpp SacFile.o string_utils.o time_utils.o
	$(CC) $(CFLAGS) $(INCBOOST) -o $@ $^ $(LIBBOOST) -lm

noise_preprocess: noise_preprocess.o SacFile.o whiten.o iir.o grbeam.o
	$(FC) $(FFLAGS) $(INCBOOST) -o $@ $^ $(LIBBOOST) -lstdc++ $(LGSL)

noise_xcorr:: noise_xcorr.o string_utils.o SacFile.o fourier.o time_utils.o
	$(CC) -static $(CFLAGS) $(INCBOOST) -o $@ $^ $(LIBBOOST) $(LFFTW) -lm

clean:
	rm -f $(PROGS) *.o core

install:
	for i in $(PROGS); do $(INSTALL) -c -m 0755 $$i $(BINDIR); done
