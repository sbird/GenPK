PRO=#-pg
#icc
# CC=icc
# CXX = icpc
# CFLAGS=-O2 -g -c -std=c99 -w1 -openmp -I/home/spb41/Lyman-alpha/GadgetReader
# CXXFLAGS= ${CFLAGS}
# LINK=${CXX} -openmp
# LFLAGS = -lfftw3f_threads -lfftw3f -lpthread -lgadread -L/home/spb41/Lyman-alpha/GadgetReader
#gcc
GREAD=/home/spb41/personal-code/GadgetReader
CC=gcc
CXX = g++
CFLAGS=-O2 -g -c -Wall -fopenmp -I${GREAD}
CXXFLAGS:= ${CFLAGS}
CFLAGS+= -std=c99
LINK=${CXX} -openmp $(PRO)
LFLAGS = -lm -lgomp -lfftw3f_threads -lfftw3f -lpthread -lgadread -L${GREAD}
objs = powerspectrum.o fieldize.o read_fieldize.o utils.o
.PHONY:all love clean test
all:gen-pk
gen-pk: gen-pk.o ${objs}
	${LINK} ${LFLAGS} $^ -o $@
powerspectrum.o: powerspectrum.c gen-pk.h
fieldize.o:fieldize.c 
read_fieldize.o: read_fieldize.cpp gen-pk.h
utils.o: utils.cpp
gen-pk.o:gen-pk.cpp gen-pk.h 
test: test.cpp ${objs}
	${LINK}-I${GREAD} ${LFLAGS} -lboost_unit_test_framework $^ -o $@
	./test
love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	rm ${objs} gen-pk.o gen-pk
