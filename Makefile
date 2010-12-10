PRO=#-pg
GREAD=/home/spb41/Lyman-alpha/GadgetReader
#icc
CC=icc
CXX = icpc
CFLAGS=-O2 -g -c -w1 -openmp -I${GREAD}
LINK=${CXX} -openmp
LFLAGS = -lfftw3f_threads -lfftw3f -lpthread -lrgad -L${GREAD} -Wl,-rpath,$(GREAD)
#gcc
# CC=gcc
# CXX = g++
# CFLAGS=-O2 -g -c -Wall -fopenmp -I${GREAD}
# LINK=${CXX} -openmp $(PRO)
# LFLAGS = -lm -lgomp -lfftw3f_threads -lfftw3f -lpthread -lgadread -L${GREAD}
CXXFLAGS:= ${CFLAGS}
CFLAGS+= -std=c99
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
btest: test.cpp ${objs}
	${LINK} -I${GREAD} ${LFLAGS} -lboost_unit_test_framework $^ -o $@
test: btest
	@./btest
love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	-rm -f ${objs} gen-pk.o gen-pk
