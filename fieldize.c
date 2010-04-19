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
/* the "extra" switch, if set to one, will assume that the output 
 * is about to be handed to an FFTW in-place routine, 
 * and set skip the last 2 places of the each row in the last dimension
 */
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#define IL 16

#ifndef M_PI
        #define M_PI 3.1415926535897932384626433832795
#endif

int fieldize(double boxsize, int dims, float *out, int total_particles, int segment_particles, float *positions,int extra)
{
	const int dims3=pow(dims,3);
	const int fdims=2*(dims/2+extra);
	/*If extra is on, we want to leave space for FFTW 
	 * to put the extra bits, so skip a couple of places.*/
	const int dims2=fdims*dims;
	const float invrho=dims3/(float)total_particles;
	const float units=dims/boxsize;
	/* This is one over density.*/
#pragma omp parallel
	{
	#pragma omp for schedule(static, 4096) nowait
	for(int index=0;index<segment_particles;index+=IL)
	{
		float dx[3],tx[3], x[3], temp[IL][8];
		int fx[3],nex[3],temp2[IL][8];
      int il=(index+IL<segment_particles ? IL : segment_particles-index);
		for(int k=0; k<il; k++)
		{
			for(int i=0; i<3; i++)
			{
				x[i]=positions[3*(index+k)+i]*units;	
				fx[i]=floor(x[i]);
				dx[i]=x[i]-fx[i];
				tx[i]=1.0-dx[i];
				nex[i]=(fx[i]+1)%dims;
            if(nex[i]<0) nex[i]+=dims;
				fx[i]%=dims;
            if(fx[i]<0) fx[i]+=dims;
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
/* 		out[dims2*fx[0] +fdims*fx[1] + fx[2]]+=invrho*tx[0]*tx[1]*tx[2]; */
/* 		out[dims2*nex[0]+fdims*fx[1] + fx[2]]+=invrho*dx[0]*tx[1]*tx[2]; */
/* 		out[dims2*fx[0] +fdims*nex[1]+ fx[2]]+=invrho*tx[0]*dx[1]*tx[2]; */
/* 		out[dims2*nex[0]+fdims*nex[1]+ fx[2]]+=invrho*dx[0]*dx[1]*tx[2]; */
/* 		out[dims2*fx[0] +fdims*fx[1] +nex[2]]+=invrho*tx[0]*tx[1]*dx[2]; */
/* 		out[dims2*nex[0]+fdims*fx[1] +nex[2]]+=invrho*dx[0]*tx[1]*dx[2]; */
/* 		out[dims2*fx[0] +fdims*nex[1]+nex[2]]+=invrho*tx[0]*dx[1]*dx[2]; */
/* 		out[dims2*nex[0]+fdims*nex[1]+nex[2]]+=invrho*dx[0]*dx[1]*dx[2]; */
		}
	}
	}
	return 0;
}

//The window function of the CiC procedure above. Need to deconvolve this for the power spectrum.
float invwindow(int kx, int ky, int kz, int n)
{
	float iwx,iwy,iwz;
	if(!kx)
		iwx=1.0;
	else
		iwx=M_PI*kx/(n*sin(M_PI*kx/(float)n));
	if(!ky)
		iwy=1.0;
	else
		iwy=M_PI*ky/(n*sin(M_PI*ky/(float)n));
	if(!kz)
		iwz=1.0;
	else
		iwz=M_PI*kz/(n*sin(M_PI*kz/(float)n));
	return pow(iwx*iwy*iwz,2);
}
