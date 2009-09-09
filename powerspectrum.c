#include <sfftw.h>
#include <sfftw_threads.h>
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
	fftw_complex *temp;
	float *powerpriv;
	int *countpriv;
	int dims2=dims*dims;
	int dims3=dims2*dims;
	//How many bins per unit interval in k?
	int binsperunit=nrbins/(floor(sqrt(3)*abs((dims+1.0)/2.0)+1));
	//Half the bin width
	float bwth=1.0/(2.0*binsperunit);
	int psindex;
	outfield=malloc(dims*dims*dims*sizeof(fftw_complex));
	temp=malloc(dims*dims*dims*sizeof(fftw_complex));
	for(int i=0; i<dims3; i++)
	{
		outfield[i].re=field[i];
		outfield[i].im=0;
	}
	if(fftw_threads_init())
	{
			  fprintf(stderr,"Error initialising fftw threads\n");
			  exit(1);
	}
/* 	pl=rfftw3d_create_plan(dims,dims,dims,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE); */
/* 	rfftwnd_threads_one_real_to_complex(4,pl, field, NULL); */
	pl=fftw3d_create_plan(dims,dims,dims,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
	fftwnd_threads_one(4,pl,outfield,NULL);
	//Now we compute the powerspectrum in each direction.
	//FFTW is unnormalised, so we need to scale by the length of the array (we do this later).
	//We could possibly dispense with some of this memory by reusing outfield, but hey.
	for(int i=0; i< nrbins; i++)
	{
		power[i]=0;
		count[i]=0;
		keffs[i]=0;
	}
	//Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
	//Use the symmetry of the real fourier transform to half the final dimension.
/* 	printf("npart=%d, cube root of npart=%d\n",npart,(int)cbrt(npart)); */
	for(int i=0; i<dims;i++)
	{
		int indx=i*dims2;
		for(int j=0; j<dims; j++)
		{
			int indy=j*dims;
			for(int k=0; k<dims; k++)
			{
				int index=indx+indy+k;
				float kk=sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2));
				psindex=floor(binsperunit*kk);
/* 				if(psindex < 2) */
/* 						  printf("%e %d %d %d\n",kk,KVAL(i),KVAL(j),KVAL(k)); */
				//Correct for shot noise.
				power[psindex]+=pow(2*M_PI,3)*(pow(outfield[index].re,2)+pow(outfield[index].im,2))/pow(dims3,2);//*pow(invwindow(KVAL(i),KVAL(j),KVAL(k),dims),2)-1.0/(double)npart;
				count[psindex]++;
			}
		}
	}
	for(int i=0; i< nrbins;i++)
	{
		
/* 	power[i]/=pow(dims3,2); */
		//bin center (k) is i+a.
		//a is bin width/2, is 0.5
		//k_eff is k+ 2a^2k/(a^2+3k^2)
		if(count[i])
		{
			float k=i*2.0*bwth;
			keffs[i]=(k+bwth)+2*pow(bwth,2)*(k+bwth)/(pow(bwth,2)+3*pow((k+bwth),2));
			power[i]/=count[i];
		}
	}
	fftwnd_destroy_plan(pl);
	free(outfield);
/* 	return 0; */
	return nrbins;

//This doesn't work.
#if 0
#pragma omp parallel private(powerpriv,countpriv)
{
	powerpriv=malloc(nrbins*sizeof(float));
	countpriv=malloc(nrbins*sizeof(int));
	if(!powerpriv || !countpriv)
	{
		fprintf(stderr,"Error allocating memory for powerspectrum!");
		exit(1);
	}
	for(int i=0; i<nrbins; i++)
	{
		powerpriv[i]=0;
		countpriv[i]=0;
	}
	#pragma omp for schedule(static)
	for(int i=0; i<dims; i++)
	{
		int indx=i*dims2;
		for(int j=0; j<dims; j++)
		{
			//k=0 case.
			int indy=j*dims;
			psindex=floor(sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)));
			powerpriv[psindex]+=(pow(outfield[indx+indy].re,2)+pow(outfield[indx+indy].im,2));
			countpriv[psindex]++;
			for(int k=0; k<dims/2; k++)
			{
				int index=indx+indy+k;
				psindex=floor(sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+pow(KVAL(k),2)));
				powerpriv[psindex]+=2*(pow(outfield[index].re,2)+pow(outfield[index].im,2));
				countpriv[psindex]+=2;
			}
			//The Nyquist frequency!
			psindex=floor(sqrt(pow(KVAL(i),2)+pow(KVAL(j),2)+dims2/4));
			powerpriv[psindex]+=(pow(outfield[indy+indx+dims/2].re,2)+pow(outfield[indx+indy+dims/2].im,2));
			countpriv[psindex]++;
		}
	}
	//Sync needs to be done carefully.
	#pragma omp critical
	{
		for(int i=0; i<nrbins;i++)
		{
			power[i]+=powerpriv[i];
			count[i]+=countpriv[i];
		}
	}
	free(powerpriv);
	free(countpriv);
	#pragma omp for schedule(static)
	//Perform normalization.
	for(int i=0; i< nrbins;i++)
	{
		power[i]/=pow(dims3,2);
		if(count[i]) power[i]/=count[i];
	}
}
	fftwnd_destroy_plan(pl);
	return nrbins;
#endif
}

