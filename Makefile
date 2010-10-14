PRO=#-pg
OF=
#-DOLD_FORMAT
#CC=gcc -O3 -g -Wall -c -std=gnu99 -fopenmp $(OF) $(PRO)
#LINK=gcc -lm -lsrfftw_threads -lsfftw_threads -lsrfftw -lsfftw -lpthread -lgomp -L/data/store/spb41/apps/fftw/lib $(PRO)
# icc; segfaults or fails to read.
CC=icc
CXX = icpc
CFLAGS=-O2 -g -c -std=c99 -w1 -openmp -I/home/spb41/Lyman-alpha/GadgetReader
CXXFLAGS= ${CFLAGS}
LINK=${CXX} -openmp $(PRO)
LFLAGS = -lsrfftw_threads -lsfftw_threads -lsrfftw -lsfftw -lpthread -lgadread -L/home/spb41/Lyman-alpha/GadgetReader
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
