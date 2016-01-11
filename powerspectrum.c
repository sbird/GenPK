/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/** \file 
 * Defines powerspectrum() wrapper around FFTW*/
extern double invwindow(int kx, int ky, int kz, int n);

/*Note we need some contiguous memory space after the actual data in field. *The real input data has size
 *dims*dims*dims
 *The output has size dims*dims*(dims/2+1) *complex* values
 * So need dims*dims*dims+2 float space.
 * Also the field needs to be stored carefully to make the 
 * extra space be in the right place. */

/**Little macro to work the storage order of the FFT.*/
#define KVAL(n) ((n)<=dims/2 ? (n) : ((n)-dims))

int powerspectrum(int dims, fftw_complex *outfield, fftw_complex *outfield2, int nrbins, double *power, int *count,double *keffs, double total_mass, double total_mass2)
{
        /*How many bins per unit (log) interval in k?*/
        const double binsperunit=(nrbins-1)/log(sqrt(3)*dims/2.0);
        /* Now we compute the powerspectrum in each direction.
         * FFTW is unnormalised, so we need to scale by the length of the array
         * (we do this later). */
        memset(power, 0, nrbins*sizeof(double));
        memset(count, 0, nrbins*sizeof(int));
        memset(keffs, 0, nrbins*sizeof(double));
        #pragma omp parallel 
        {
                double powerpriv[nrbins], keffspriv[nrbins];
                int countpriv[nrbins];
                memset(powerpriv, 0, nrbins*sizeof(double));
                memset(keffspriv, 0, nrbins*sizeof(double));
                memset(countpriv, 0, nrbins*sizeof(int));
                /* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
                 * Use the symmetry of the real fourier transform to half the final dimension.*/
                #pragma omp for nowait
                for(int i=0; i<dims;i++){
                        int indx=i*dims*(dims/2+1);
                        for(int j=0; j<dims; j++){
                                int indy=j*(dims/2+1);
                                /* The k=0 and N/2 mode need special treatment here, 
                                 * as they alone are not doubled.*/
                                /*Do k=0 mode.*/
                                int index=indx+indy;
                                double kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2));
                                //We don't want the 0,0,0 mode as that is just the mean of the field.
                                if (kk > 0) {
                                    int psindex=floor(binsperunit*log(kk));
                                    assert(psindex < nrbins);
                                    powerpriv[psindex] += (outfield[index][0]*outfield2[index][0]+outfield[index][1]*outfield2[index][1])*pow(invwindow(KVAL(i),KVAL(j),0,dims),2);
                                    keffspriv[psindex]+=kk;
                                    countpriv[psindex]++;
                                }
                                /*Now do the k=N/2 mode*/
                                index=indx+indy+dims/2;
                                kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(dims/2),2));
                                int psindex=floor(binsperunit*log(kk));
                                assert(psindex < nrbins);
                                powerpriv[psindex] += (outfield[index][0]*outfield2[index][0]+outfield[index][1]*outfield2[index][1])*pow(invwindow(KVAL(i),KVAL(j),KVAL(dims/2),dims),2);
                                keffspriv[psindex]+=kk;
                                countpriv[psindex]++;
                                /*Now do the rest. Because of the symmetry, each mode counts twice.*/
                                for(int k=1; k<dims/2; k++){
                                        index=indx+indy+k;
                                        kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2));
                                        int psindex=floor(binsperunit*log(kk));
                                        assert(psindex < nrbins);
                                        powerpriv[psindex]+=2*(outfield[index][0]*outfield2[index][0]+outfield[index][1]*outfield2[index][1])*pow(invwindow(KVAL(i),KVAL(j),KVAL(k),dims),2);
                                        countpriv[psindex]+=2;
                                        keffspriv[psindex]+=2*kk;
                                }
                        }
                }
                //Can't do reductions on arrays yet.
                #pragma omp critical
                {
                        for(int i=0; i< nrbins;i++){
                                power[i]+=powerpriv[i];
                                count[i]+=countpriv[i];
                                keffs[i]+=keffspriv[i];
                        }
                }
        }
        for(int i=0; i< nrbins;i++){
                if(count[i]){
                        power[i]/=total_mass*total_mass2;
                        power[i]/=count[i];
                        keffs[i]/=count[i];
                }
        }
        return 0;
}

