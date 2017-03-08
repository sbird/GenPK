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

/* Fieldize. positions should be an array of size 3*segment_particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
#include <cmath>
#include <stdint.h>
#include "gen-pk.h"
#include <cassert>

/** Number of particles to keep in the thread-local buffer*/
#define IL 64

#ifndef M_PI
        /** C99 does not define PI!*/
        #define M_PI 3.1415926535897932384626433832795
#endif

/** \file
 * Defines fieldize() 
 * Fieldize: this does cloud-in-cell interpolation, assuming periodic boundary conditions. 
 * @param boxsize The physical size of the grid. This is only used so to convert the returned data to physical units
 * @param dims The number of grid points to use
 * @param out Pointer to the grid to interpolate onto. 
 *          Should be an array of size dims^3 (or slightly larger if extra is set, see below.
 * @param segment_particles The number of particles to interpolate on this call. The length of the positions array.
 * @param positions Array of particle positions. Assumed to have dimensions 3*segment_particles
 * @param masses Array containing particle masses, if variable. Null if particle masses are constant.
 * @param mass Particle mass if constant for all of this type.
 * @param extra This switch should be 0 or 1. If set to one, will assume that the output
 * is about to be handed to an FFTW in-place routine, 
 * and skip the last 2 places of each row in the last dimension (of out)
 */
int fieldize(double boxsize, int dims, GENFLOAT *out, int64_t segment_particles, GENPK_FLOAT_TYPE *positions, GENPK_FLOAT_TYPE * masses, double mass, int extra)
{
    const size_t fdims=2*(dims/2+extra);
    /*If extra is on, we want to leave space for FFTW
     * to put the extra bits, so skip a couple of places.*/
    const size_t dims2=fdims*dims;
    const GENFLOAT units=dims/boxsize;
    #pragma omp parallel for
    for(int64_t index=0;index<segment_particles;index+=IL) {
        GENFLOAT temp[IL][8];
        size_t temp2[IL][8];
        const int il=(index+IL<segment_particles ? IL : segment_particles-index);
        for(int k=0; k<il; k++)
        {
            GENFLOAT dx[3],tx[3], x[3];
            int fx[3],nex[3];
            /* This is one over density.*/
            const GENFLOAT invrho = (masses ? masses[index+k] : mass);
            for(int i=0; i<3; i++)
            {
                x[i]=positions[3*(index+k)+i]*units;
                fx[i]=floor(x[i]);
                dx[i]=x[i]-fx[i];
                tx[i]=1.0-dx[i];
                nex[i]=(fx[i]+1)%dims;
                if(nex[i]<0) 
                    nex[i]+=dims;
                fx[i]%=dims;
                if(fx[i]<0)
                    fx[i]+=dims;
            }
            temp[k][0]=invrho*tx[0]*tx[1]*tx[2];
            temp[k][1]=invrho*dx[0]*tx[1]*tx[2];
            temp[k][2]=invrho*tx[0]*dx[1]*tx[2];
            temp[k][3]=invrho*dx[0]*dx[1]*tx[2];
            temp[k][4]=invrho*tx[0]*tx[1]*dx[2];
            temp[k][5]=invrho*dx[0]*tx[1]*dx[2];
            temp[k][6]=invrho*tx[0]*dx[1]*dx[2];
            temp[k][7]=invrho*dx[0]*dx[1]*dx[2];
            temp2[k][0]=dims2*fx[0] +fdims*fx[1] + fx[2];
            temp2[k][1]=dims2*nex[0]+fdims*fx[1] + fx[2];
            temp2[k][2]=dims2*fx[0] +fdims*nex[1]+ fx[2];
            temp2[k][3]=dims2*nex[0]+fdims*nex[1]+ fx[2];
            temp2[k][4]=dims2*fx[0] +fdims*fx[1] +nex[2];
            temp2[k][5]=dims2*nex[0]+fdims*fx[1] +nex[2];
            temp2[k][6]=dims2*fx[0] +fdims*nex[1]+nex[2];
            temp2[k][7]=dims2*nex[0]+fdims*nex[1]+nex[2];
            for(int c=0; c<8; c++)
                assert(temp2[k][c] >= 0 && temp2[k][c] < dims2*dims);
        }
        /*The store operation may only be done by one thread at a time,
        *to ensure synchronisation.*/
        #pragma omp critical
        {
          for(int k=0; k<il; k++)
          {
            out[temp2[k][0]]+=temp[k][0];
            out[temp2[k][1]]+=temp[k][1];
            out[temp2[k][2]]+=temp[k][2];
            out[temp2[k][3]]+=temp[k][3];
            out[temp2[k][4]]+=temp[k][4];
            out[temp2[k][5]]+=temp[k][5];
            out[temp2[k][6]]+=temp[k][6];
            out[temp2[k][7]]+=temp[k][7];
          }
        }
    }
    return 0;
}

//Helper function for 1D window function.
inline GENFLOAT onedinvwindow(int kx, int n)
{
    //Return \pi x /(n sin(\pi x n)) unless x = 0, in which case return 1.
    return kx ? M_PI*kx/(n*sin(M_PI*kx/(float)n)) : 1.0;
}

//The window function of the CiC procedure above. Need to deconvolve this for the power spectrum.
//Only has an effect for k > Nyquist/4.
GENFLOAT invwindow(int kx, int ky, int kz, int n)
{
    if(n == 0)
        return 0;
    float iwx = onedinvwindow(kx, n);
    float iwy = onedinvwindow(ky, n);
    float iwz = onedinvwindow(kz, n);
    return pow(iwx*iwy*iwz,2);
}
