/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
//The overhead involved in openmp makes the parallel version about 
//5 times slower than the serial. Awesome.
//#include <omp.h>
#include <math.h>

int fieldize(float boxsize, int dims, float *out, int particles, float *positions)
{
	int dims3=pow(dims,3);
	int dims2=pow(dims,2);
	float invrho=dims3/particles;
	int fx[3],nex[3],i,index;
	float dx[3],tx[3], temp[8];
	float x[3];
	float units=dims/boxsize;
	/* This is one over density.*/
/* 	#pragma omp parallel private(i,index,temp) */
/* 	{ */
/* 		#pragma omp for schedule(static) */
		for(int i=0; i<dims*dims*dims; i++)
			out[i]=-1;
/* 		#pragma omp for schedule(static) */
		for(int index=0;index<particles;index++)
		{
			for(i=0; i<3; i++)
			{
				x[i]=positions[3*index+i]*units;
				fx[i]=floor(x[i]);
				dx[i]=x[i]-fx[i];
				tx[i]=1.0-dx[i];
				nex[i]=(fx[i]+1)%dims;
				fx[i]%=dims;
			}
/*			temp[0]=invrho*tx[0]*tx[1]*tx[2];
			temp[1]=invrho*dx[0]*tx[1]*tx[2];
			temp[2]=invrho*tx[0]*dx[1]*tx[2];
			temp[3]=invrho*dx[0]*dx[1]*tx[2];
			temp[4]=invrho*tx[0]*tx[1]*dx[2];
			temp[5]=invrho*dx[0]*tx[1]*dx[2];
			temp[6]=invrho*tx[0]*dx[1]*dx[2];
			temp[7]=invrho*dx[0]*dx[1]*dx[2];*/
/* 			#pragma omp critical */
/* 			{  */
		//The store operation may only be done by one thread at a time, 
		//to ensure synchronisation.
			out[dims2*fx[0] +dims*fx[1] +fx[2]]	+=invrho*tx[0]*tx[1]*tx[2];
			out[dims2*nex[0]+dims*fx[1] +fx[2]]	+=invrho*dx[0]*tx[1]*tx[2];
			out[dims2*fx[0] +dims*nex[1]+fx[2]]	+=invrho*tx[0]*dx[1]*tx[2];
			out[dims2*nex[0]+dims*nex[1]+fx[2]]	+=invrho*dx[0]*dx[1]*tx[2];
			out[dims2*fx[0] +dims*fx[1] +nex[2]]+=invrho*tx[0]*tx[1]*dx[2];
			out[dims2*nex[0]+dims*fx[1] +nex[2]]+=invrho*dx[0]*tx[1]*dx[2];
			out[dims2*fx[0] +dims*nex[1]+nex[2]]+=invrho*tx[0]*dx[1]*dx[2];
			out[dims2*nex[0]+dims*nex[1]+nex[2]]+=invrho*dx[0]*dx[1]*dx[2];
/* 			} */
		}
/* 	} */
	return 0;
}
