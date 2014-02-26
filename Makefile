
#Change this to where you installed GadgetReader
GREAD=${CURDIR}/../GadgetReader


ifeq ($(CC),cc)
  ICC:=$(shell which icc --tty-only 2>&1)
  #Can we find icc?
  ifeq (/icc,$(findstring /icc,${ICC}))
     CC = icc -vec_report0
     CXX = icpc
  else
     GCC:=$(shell which gcc --tty-only 2>&1)
     #Can we find gcc?
     ifeq (/gcc,$(findstring /gcc,${GCC}))
        CC = gcc
        CXX = g++
     endif
  endif
endif

LFLAGS += -lfftw3f_threads -lfftw3f -lpthread -lrgad -L${GREAD} -Wl,-rpath,$(GREAD) -lhdf5 -lhdf5_hl
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp -I${GREAD}
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -g -c -Wall -fopenmp -I${GREAD}
  LINK +=${CXX} -openmp $(PRO)
  LFLAGS += -lm -lgomp 
endif
PRO=#-pg
#gcc
PPFLAGS:=$(CFLAGS)
CXXFLAGS+= $(PPFLAGS)
CFLAGS+= -std=c99
objs = powerspectrum.o fieldize.o read_fieldize.o utils.o
.PHONY:all love clean test dist

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

dist: Makefile README $(head) Doxyfile gen-pk.cpp  read_fieldize.cpp  test.cpp  utils.cpp gen-pk.h fieldize.c powerspectrum.c test_g2_snap.0 test_g2_snap.1
	tar -czf genpk.tar.gz $^

doc: Doxyfile gen-pk.h
	doxygen $<

love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	-rm -f ${objs} gen-pk.o gen-pk
