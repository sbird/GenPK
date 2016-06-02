
#Change this to where you installed GadgetReader
GREAD=${CURDIR}/GadgetReader

#If Gadget was compiled with double precision output, you should define this flag
#to read it correctly
#OPT = -DDOUBLE_PRECISION

LFLAGS += -lfftw3_threads -lfftw3 -lpthread -lrgad -L${GREAD} -Wl,-rpath,$(GREAD),--no-add-needed,--as-needed -lhdf5 -lhdf5_hl
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp $(OPT) -I${GREAD}
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -ffast-math -g -c -Wall -fopenmp $(OPT) -I${GREAD}
  LINK +=${CXX} -openmp $(PRO)
  LFLAGS += -lm -lgomp 
endif
PRO=#-pg
#gcc
PPFLAGS:=$(CFLAGS)
CXXFLAGS+= $(PPFLAGS)
objs = powerspectrum.o fieldize.o read_fieldize.o utils.o
.PHONY:all love clean test dist

all: librgad.so gen-pk

gen-pk: gen-pk.o ${objs}
	${LINK} $^ ${LFLAGS} -o $@

librgad.so:
	cd $(GREAD); VPATH=$(GREAD) make $@

powerspectrum.o: powerspectrum.c
	$(CC) -std=gnu99 $(CFLAGS) $^
%.o: %.cpp gen-pk.h

btest: test.cpp ${objs}
	${LINK} $(OPT) -I${GREAD} $^ ${LFLAGS} -lboost_unit_test_framework -o $@

test: btest librgad.so
	./$<

dist: Makefile README $(head) Doxyfile gen-pk.cpp  read_fieldize.cpp  test.cpp  utils.cpp gen-pk.h fieldize.cpp powerspectrum.c test_g2_snap.0 test_g2_snap.1
	tar -czf genpk.tar.gz $^

doc: Doxyfile gen-pk.h
	doxygen $<

clean:
	-rm -f ${objs} gen-pk.o gen-pk
