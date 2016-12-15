
#Change this to where you installed GadgetReader
GREAD=${CURDIR}/GadgetReader

#If Gadget was compiled with double precision output, you should define this flag
#to read it correctly
#OPT = -DDOUBLE_PRECISION

#Check for a pkgconfig; if one exists we are probably debian.
ifeq ($(shell pkg-config --exists hdf5-serial && echo 1),1)
	HDF_LINK = $(shell pkg-config --libs hdf5-serial) -lhdf5_hl
	HDF_INC = $(shell pkg-config --cflags hdf5-serial)
else
	HDF_LINK = -lhdf5 -lhdf5_hl
endif

LFLAGS += -lfftw3_threads -lfftw3 -lpthread -lrgad -L${GREAD} -Wl,-rpath,$(GREAD),--no-add-needed,--as-needed $(HDF_LINK) -Lbigfile/src -lbigfile
#Are we using gcc or icc?
ifeq (icc,$(findstring icc,${CC}))
  CFLAGS +=-O2 -g -c -w1 -openmp $(OPT) -I${GREAD}
  LINK +=${CXX} -openmp
else
  CFLAGS +=-O2 -ffast-math -g -c -Wall -fopenmp $(OPT) -I${GREAD} $(HDF_INC)
  LINK +=${CXX} -openmp $(PRO)
  LFLAGS += -lm -lgomp 
endif
PRO=#-pg
#gcc
PPFLAGS:=$(CFLAGS)
CXXFLAGS+= $(PPFLAGS)
objs = powerspectrum.o fieldize.o read_fieldize.o utils.o read_fieldize_bigfile.o
.PHONY:all love clean test dist

all: librgad.so gen-pk

gen-pk: gen-pk.o ${objs}
	${LINK} $^ ${LFLAGS} -o $@

$(GREAD)/Makefile:
	git submodule init
	git submodule update

librgad.so: $(GREAD)/Makefile
	cd $(GREAD); VPATH=$(GREAD) make $@

read_fieldize_bigfile.o: bigfile/src/libbigfile.a

bigfile/src/libbigfile.a:
	cd bigfile/src; VPATH=bigfile/src MPICC=$(CC) make libbigfile.a

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
