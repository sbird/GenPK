/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
/* the "extra" switch, if set to one, will assume that the output 
 * is about to be handed to an FFTW in-place routine, 
 * and set skip the last 2 places of the each row in the last dimension
 */
#include <math.h>
#include <stdio.h>
#define IL 16

int fieldize(double boxsize, int dims, float *out, int particles, float *positions,int extra)
{
	const int dims3=pow(dims,3);
	const int fdims=2*(dims/2+extra);
	/*If extra is on, we want to leave space for FFTW 
	 * to put the extra bits, so skip a couple of places.*/
	const int dims2=fdims*dims;
	const float invrho=dims3/(float)particles;
	const float units=dims/boxsize;
	/* This is one over density.*/
#pragma omp parallel
	{
	#pragma omp for schedule(static, 1024)
	for(int i=0; i<dims2; i+=fdims)
		for(int j=0; j<dims;j++)
				out[i+j]=-1;
	#pragma omp for schedule(static, 32768) nowait
	for(int index=0;index<particles;index+=IL)
	{
		float dx[3],tx[3], x[3], temp[IL][8];
		int fx[3],nex[3],temp2[IL][8];
		for(int k=0; k<IL; k++)
		{
			for(int i=0; i<3; i++)
			{
				x[i]=positions[3*(index+k)+i]*units;	
				fx[i]=floor(x[i]);
				dx[i]=x[i]-fx[i];
				tx[i]=1.0-dx[i];
				nex[i]=(fx[i]+1)%dims;
				fx[i]%=dims;
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
		for(int k=0; k<IL; k++)
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
