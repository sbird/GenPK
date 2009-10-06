PRO=-pg
OF=
#-DOLD_FORMAT
CC=gcc -O3  -c -std=gnu99 $(OF) -fopenmp $(PRO)
LINK=gcc -lm -lsrfftw_threads -lsfftw_threads -lsrfftw -lsfftw -lpthread -lgomp -L/data/store/spb41/apps/fftw/lib $(PRO)
# icc; segfaults or fails to read.
#CC=icc -O2 -g -c -c99
#LINK=icc -lm -lsrfftw_threads -lsfftw_threads -lsrfftw -lsfftw -lpthread -L/data/store/spb41/apps/fftw/lib $(PRO)
.PHONY:all love clean
all:gen-pk
gen-pk:powerspectrum.o fieldize.o readgadget.o gen-pk.c gen-pk.h
	${CC} gen-pk.c
	${LINK} fieldize.o powerspectrum.o readgadget.o gen-pk.o -o gen-pk
powerspectrum.o: powerspectrum.c
	${CC} powerspectrum.c
fieldize.o:fieldize.c
	${CC} fieldize.c
readgadget.o: readgadget.c readgadget.h
	${CC} readgadget.c
love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	rm *.o
