PRO=#-pg
#icc
CC=icc
CXX = icpc
CFLAGS=-O2 -g -c -std=c99 -w1 -openmp -I/home/spb41/Lyman-alpha/GadgetReader
CXXFLAGS= ${CFLAGS}
LINK=${CXX} -openmp
LFLAGS = -lfftw3f_threads -lfftw3f -lpthread -lgadread -L/home/spb41/Lyman-alpha/GadgetReader
#gcc
# CC=gcc
# CXX = g++
# CFLAGS=-O0 -g -c -Wall -fopenmp -I/home/spb41/Lyman-alpha/GadgetReader
# CFLAGS+= -std=c99
# LINK=${CXX} -openmp $(PRO)
# LFLAGS = -lm -lgomp -lfftw3f_threads -lfftw3f -lpthread -lgadread -L/home/spb41/Lyman-alpha/GadgetReader
objs = gen-pk.o powerspectrum.o fieldize.o 
.PHONY:all love clean
all:gen-pk
gen-pk: ${objs}
	${LINK} ${LFLAGS} $^ -o $@
powerspectrum.o: powerspectrum.c 
fieldize.o:fieldize.c 
gen-pk.o:gen-pk.cpp gen-pk.h 
love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	rm ${objs} gen-pk
