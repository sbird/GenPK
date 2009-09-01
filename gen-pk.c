#include <math.h>
#include <stdio.h>
#include "readgadget.h"
#include <sfftw.h>
/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
int fieldize(float boxsize, int dims, float *out, int particles, float *positions);

/* The window function associated with the above.*/
float invwindow(int kx, int ky, int kz, int n);
//Note we will need some contiguous memory space after the actual data in field.
// The real input data has size
// dims*dims*dims
// The output has size dims*dims*(dims/2+1) *complex* values
// So need 2*(dims*dims*dims+2) float space.
// Need at least floor(sqrt(3)*abs((dims+1.0)/2.0)+1) values in power and count.
int powerspectrum(int dims, fftw_real *field, float *power, float *count, float *keffs, int particles);

#define FIELD_DIMS 128

int main(char argc, char* argv[]){
  int npart[5],nrbins=2*floor(sqrt(3)*abs((FIELD_DIMS+1.0)/2.0)+1);
  double massarr[5], time, redshift;
  float boxsize;
  float *pos, *power, *count, *keffs;
  float *field;
  if(argc<3)
  {
			 printf("Need filename and boxsize\n");
			 exit(0);
  }
  boxsize=atof(argv[2]);
  FILE *fd=fopen(argv[1],"r");
  if(!fd)
  {
		fprintf(stderr,"Error opening file for reading!\n");
		exit(1);
  }
  if(!read_gadget_head(npart, massarr, &time, &redshift, fd))
  {
		fprintf(stderr,"Error reading file header!\n");
		exit(1);
  }
  pos=malloc(3*(npart[1]+1)*sizeof(float));
  field=malloc((FIELD_DIMS*FIELD_DIMS*FIELD_DIMS+2)*sizeof(float));
  if(!pos | !field)
  {
		fprintf(stderr,"Error allocating particle memory\n");
		exit(1);
  }
  if(!read_gadget_float3(pos, "POS ",fd))
  {
		fprintf(stderr, "Error reading particle data\n");
		exit(1);
  }
  //By now we should have the data.
  fclose(fd);
  fieldize(boxsize,FIELD_DIMS,field,npart[1],pos);
  free(pos);
  power=malloc(nrbins*sizeof(float));
  count=malloc(nrbins*sizeof(float));
  keffs=malloc(nrbins*sizeof(float));
  if(!power | !count)
  {
		fprintf(stderr,"Error allocating memory for power spectrum.\n");
		exit(1);
  }
  nrbins=powerspectrum(FIELD_DIMS,(fftw_real *)field,power,count,keffs,npart[1]);
  free(field);
 for(int i=0;i<nrbins;i++)
 {
	printf("%e\t%e\t%e\t%e\n",i+0.5,keffs[i],power[i],count[i]);
 }
	free(power);
	free(count);
  return 0;
}
