/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
//The overhead involved in openmp makes the parallel version about 
//5 times slower than the serial. Awesome.
#include <math.h>
#include <stdio.h>

int fieldize(double boxsize, int dims, float *out, int particles, float *positions)
{
	const int dims3=pow(dims,3);
	const int dims2=pow(dims,2);
	const float invrho=dims3/(float)particles;
	const float units=dims/boxsize;
	/* This is one over density.*/
#pragma omp parallel
	{
	#pragma omp for schedule(static)
	for(int i=0; i<dims3; i++)
		out[i]=-1;
	#pragma omp for schedule(static) nowait
	for(int index=0;index<particles;index++)
	{
		float dx[3],tx[3], x[3];
		int fx[3],nex[3];
		for(int i=0; i<3; i++)
		{
			x[i]=positions[3*index+i]*units;	
			fx[i]=floor(x[i]);
			dx[i]=x[i]-fx[i];
			tx[i]=1.0-dx[i];
			nex[i]=(fx[i]+1)%dims;
			fx[i]%=dims;
		}
		//The store operation may only be done by one thread at a time, 
		//to ensure synchronisation.
		#pragma omp critical 
		{
		out[dims2*fx[0] +dims*fx[1] +fx[2]]	+=invrho*tx[0]*tx[1]*tx[2];
		out[dims2*nex[0]+dims*fx[1] +fx[2]]	+=invrho*dx[0]*tx[1]*tx[2];
		out[dims2*fx[0] +dims*nex[1]+fx[2]]	+=invrho*tx[0]*dx[1]*tx[2];
		out[dims2*nex[0]+dims*nex[1]+fx[2]]	+=invrho*dx[0]*dx[1]*tx[2];
		out[dims2*fx[0] +dims*fx[1] +nex[2]]+=invrho*tx[0]*tx[1]*dx[2];
		out[dims2*nex[0]+dims*fx[1] +nex[2]]+=invrho*dx[0]*tx[1]*dx[2];
		out[dims2*fx[0] +dims*nex[1]+nex[2]]+=invrho*tx[0]*dx[1]*dx[2];
		out[dims2*nex[0]+dims*nex[1]+nex[2]]+=invrho*dx[0]*dx[1]*dx[2];
		}
	}
	}
	return 0;
}

//The window function of the CiC prodcedure above. Need to deconvolve this for the power spectrum.
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
