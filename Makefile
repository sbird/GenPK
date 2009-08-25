CC=gcc -O2 -g -c -std=c99
LINK=gcc -lm -lsrfftw_threads -lsfftw_threads -lsrfftw -lsfftw -lpthread -L/data/store/spb41/apps/fftw/lib 

.PHONY:all love clean
all:gen-pk

gen-pk:powerspectrum.o fieldize.o readgadget.o gen-pk.c
	${CC} gen-pk.c
	${LINK} fieldize.o powerspectrum.o readgadget.o gen-pk.o -o gen-pk
powerspectrum.o: powerspectrum.c
	${CC} powerspectrum.c
fieldize.o:fieldize.c
	${CC} fieldize.c
readgadget.o: readgadget.c
	${CC} readgadget.c
love:
	@echo "Not war?" ; sleep 3
	@echo "Look, I'm not equipped for that, okay?" ; sleep 2
	@echo "Contact your hardware vendor for appropriate mods."

clean:
	rm *.o
