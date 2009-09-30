#include "gen-pk.h"
#include <omp.h>

/*Note we need some contiguous memory space after the actual data in field. *The real input data has size
 *dims*dims*dims
 *The output has size dims*dims*(dims/2+1) *complex* values
 * So need dims*dims*dims+2 float space.
 * Also the field needs to be stored carefully to make the 
 * extra space be in the right place. */

/*Little macro to work the storage order of the FFT.*/
#define KVAL(n) ((n)<=dims/2 ? (n) : ((n)-dims))

int powerspectrum(int dims, float *field, int nrbins, float *power, float *count,float *keffs)
{
	fftwnd_plan pl;
	fftw_complex *outfield;
	const int dims2=dims*dims;
	const int dims3=dims2*dims;
	/*How many bins per unit interval in k?*/
	const int binsperunit=nrbins/(floor(sqrt(3)*abs((dims+1.0)/2.0)+1));
	/*Half the bin width*/
	const float bwth=1.0/(2.0*binsperunit);
	if(sizeof(fftw_real) != sizeof(float))
	{
		fprintf(stderr, "sizeof fftw_real:%d fftw_complex: %d, float: %d\n",sizeof(fftw_real), sizeof(fftw_complex), sizeof(float));
		fprintf(stderr, "fftw_real is not a float. Perhaps you linked the wrong library?\n");
		exit(1);
	}
	outfield=(fftw_complex *) &field[0];
	if(!outfield){
			  fprintf(stderr, "Error allocating memory for outfield!\n");
			  exit(1);
	}
	if(fftw_threads_init())
	{
			  fprintf(stderr,"Error initialising fftw threads\n");
			  exit(1);
	}
	pl=rfftw3d_create_plan(dims,dims,dims,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
	rfftwnd_threads_one_real_to_complex(omp_get_num_procs(),pl, &field[0], outfield);
	/* Now we compute the powerspectrum in each direction.
	 * FFTW is unnormalised, so we need to scale by the length of the array
	 * (we do this later). */
	for(int i=0; i< nrbins/2; i++)
	{
		/* bin center (k) is i+a.
		 * a is bin width/2, is 0.5
		 * k_eff is k+ 2a^2k/(a^2+3k^2) */
		float k=i*2.0*bwth;
		keffs[i]=(k+bwth)+2*pow(bwth,2)*(k+bwth)/(pow(bwth,2)+3*pow((k+bwth),2));
		power[i]=0;
		count[i]=0;
	}
	/*After this point, the number of modes is decreasing.*/
	for(int i=nrbins/2; i< nrbins; i++)
	{
		/* bin center (k) is i+a.
		 * a is bin width/2, is 0.5
		 * k_eff is k+ 2a^2k/(a^2+3k^2) */
		float k=i*2.0*bwth;
		keffs[i]=(k+bwth)-2*pow(bwth,2)*(k+bwth)/(pow(bwth,2)+3*pow((k+bwth),2));
		power[i]=0;
		count[i]=0;
	}
	#pragma omp parallel 
	{
		float powerpriv[nrbins];
		int countpriv[nrbins];
		for(int i=0; i< nrbins; i++)
		{
			powerpriv[i]=0;
			countpriv[i]=0;
		}
		/* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
		 * Use the symmetry of the real fourier transform to half the final dimension.*/
		#pragma omp for schedule(static, 128) nowait
		for(int i=0; i<dims;i++)
		{
			int indx=i*dims*(dims/2+1);
			for(int j=0; j<dims; j++)
			{
				int indy=j*(dims/2+1);
				/* The k=0 and N/2 mode need special treatment here, 
				 * as they alone are not doubled.*/
				/*Do k=0 mode.*/
				int index=indx+indy;
				float kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2));
				int psindex=floor(binsperunit*kk);
				powerpriv[psindex]+=(pow(outfield[index].re,2)+pow(outfield[index].im,2))*pow(invwindow(KVAL(i),KVAL(j),0,dims),2);
				countpriv[psindex]++;
				/*Now do the k=N/2 mode*/
				index=indx+indy+dims/2;
				kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(dims/2),2));
				psindex=floor(binsperunit*kk);
				powerpriv[psindex]+=(pow(outfield[index].re,2)+pow(outfield[index].im,2))*pow(invwindow(KVAL(i),KVAL(j),KVAL(dims/2),dims),2);
				countpriv[psindex]++;
				/*Now do the rest. Because of the symmetry, each mode counts twice.*/
				for(int k=1; k<dims/2; k++)
				{
					index=indx+indy+k;
					kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2));
					psindex=floor(binsperunit*kk);
					/* Correct for shot noise and window function in IDL. 
					 * See my notes for the reason why.*/
					powerpriv[psindex]+=2*(pow(outfield[index].re,2)+pow(outfield[index].im,2))*pow(invwindow(KVAL(i),KVAL(j),KVAL(k),dims),2);
					countpriv[psindex]+=2;
				}
			}
		}
		#pragma omp critical
		{
			for(int i=0; i< nrbins;i++)
			{
				power[i]+=powerpriv[i];
				count[i]+=countpriv[i];
			}
		}
	}
	for(int i=0; i< nrbins;i++)
	{
		if(count[i])
		{
			/* I do the division twice to avoid any overflow.*/
			power[i]/=dims3;
			power[i]/=dims3;
			power[i]/=count[i];
		}
	}
	fftwnd_destroy_plan(pl);
/* 	return 0; */
	return nrbins;
}

