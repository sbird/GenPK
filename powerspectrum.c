#include <srfftw.h>
#include <srfftw_threads.h>
#include <math.h>
//Note we will need some contiguous memory space after the actual data in field.
// The real input data has size
// dims*dims*dims
// The output has size dims*dims*(dims/2+1) *complex* values
// So need dims*dims*dims+2 float space.

//Little macro to work the storage order of the FFT.
#define KVAL(n) ((n)<=dims/2 ? (n) : ((n)-dims))


extern float invwindow(int kx, int ky, int kz, int n);

int powerspectrum(int dims, fftw_real *field, int nrbins, float *power, float *count,float *keffs,int npart)
{
	fftwnd_plan pl;
	fftw_complex *outfield;
	float *powerpriv;
	int *countpriv;
	int dims2=dims*dims;
	int dims3=dims2*dims;
	/*How many bins per unit interval in k?*/
	int binsperunit=nrbins/(floor(sqrt(3)*abs((dims+1.0)/2.0)+1));
	/*Half the bin width*/
	float bwth=1.0/(2.0*binsperunit);
	int psindex;
	int totalpts=dims*dims*(dims/2+1)*sizeof(fftw_real);
	if(sizeof(fftw_real) != sizeof(float))
	{
		fprintf(stderr, "sizeof fftw_real:%d fftw_complex: %d, float: %d\n",sizeof(fftw_real), sizeof(fftw_complex), sizeof(float));
		fprintf(stderr, "fftw_real is not a float. Perhaps you linked the wrong library?\n");
		exit(1);
	}
	//Need to dispense with this memory by, eg, re-using field.
	outfield=malloc(2*dims*dims*(dims/2+1)*sizeof(fftw_real));
	if(!outfield){
			  fprintf(stderr, "Error allocating memory for outfield!\n");
			  exit(1);
	}
	if(fftw_threads_init())
	{
			  fprintf(stderr,"Error initialising fftw threads\n");
			  exit(1);
	}
	pl=rfftw3d_create_plan(dims,dims,dims,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	rfftwnd_threads_one_real_to_complex(4,pl, field, outfield);
	/* Now we compute the powerspectrum in each direction.
	 * FFTW is unnormalised, so we need to scale by the length of the array
	 * (we do this later). */
	/*We could possibly dispense with some of this memory by reusing outfield, but hey.*/
	for(int i=0; i< nrbins; i++)
	{
		power[i]=0;
		count[i]=0;
		keffs[i]=0;
	}
	/* Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
	 * Use the symmetry of the real fourier transform to half the final dimension.*/
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
			psindex=floor(binsperunit*kk);
			power[psindex]+=pow(outfield[index].re,2)+pow(outfield[index].im,2);
			count[psindex]++;
			/*Now do the k=N/2 mode*/
			index=indx+indy+dims/2;
			kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(dims/2),2));
			psindex=floor(binsperunit*kk);
			power[psindex]+=pow(outfield[index].re,2)+pow(outfield[index].im,2);
			count[psindex]++;
			/*Now do the rest. Because of the symmetry, each mode counts twice.*/
			for(int k=1; k<dims/2; k++)
			{
				index=indx+indy+k;
				kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2));
				psindex=floor(binsperunit*kk);
				/* Correct for shot noise and window function in IDL. 
				 * See my notes for the reason why.*/
				power[psindex]+=2*(pow(outfield[index].re,2)+pow(outfield[index].im,2));
				count[psindex]+=2;
			}

		}
	}
	for(int i=0; i< nrbins;i++)
	{
		
		/* bin center (k) is i+a.
		 * a is bin width/2, is 0.5
		 * k_eff is k+ 2a^2k/(a^2+3k^2) */
		if(count[i])
		{
			float k=i*2.0*bwth;
			keffs[i]=(k+bwth)+2*pow(bwth,2)*(k+bwth)/(pow(bwth,2)+3*pow((k+bwth),2));
			/* I do the division twice to avoid any overflow.*/
			power[i]/=dims3;
			power[i]/=dims3;
			power[i]/=count[i];
		}
	}
	fftwnd_destroy_plan(pl);
	free(outfield);
/* 	return 0; */
	return nrbins;
}

