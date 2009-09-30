#include "gen-pk.h"
/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]
 * the "extra" switch, if set to one, will assume that the output 
 * is about to be handed to an FFTW in-place routine, 
 * and set skip the last 2 places of the each row in the last dimension
 */

/* In practice this means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 1024

int main(char argc, char* argv[]){
  int field_dims;
  int old =0;
  int npart[5],nrbins;
  double massarr[5], time, redshift;
  double boxsize;
  float *pos, *power, *count, *keffs;
  float *field;
  if(argc<2)
  {
			 printf("Need filename\n");
			 exit(0);
  }
  FILE *fd=fopen(argv[1],"r");
  if(argc>2)
  {
      old=atoi(argv[2]);
  }
  if(!fd)
  {
		fprintf(stderr,"Error opening file for reading!\n");
		exit(1);
  }
  if(!read_gadget_head(npart, massarr, &time, &redshift,&boxsize, fd, old))
  {
		fprintf(stderr,"Error reading file header!\n");
		exit(1);
  }
  fprintf(stderr, "Boxsize=%g, npart=%g, redshift=%g\n",boxsize,cbrt(npart[1]),redshift);
  field_dims=(2*cbrt(npart[1]) < FIELD_DIMS ? 2*cbrt(npart[1]) : FIELD_DIMS);
  nrbins=floor(sqrt(3)*abs((field_dims+1.0)/2.0)+1);
  pos=malloc(3*(npart[1]+1)*sizeof(float));
  /* Allocating a bit more memory allows us to do in-place transforms.*/
  field=malloc(2*field_dims*field_dims*(field_dims/2+1)*sizeof(float));
  if(!pos || !field)
  {
		fprintf(stderr,"Error allocating particle memory\n");
		exit(1);
  }
  if(!read_gadget_float3(pos, "POS ",fd, old))
  {
		fprintf(stderr, "Error reading particle data\n");
		exit(1);
  }
  //By now we should have the data.
  fclose(fd);
  fieldize(boxsize,field_dims,field,npart[1],pos, 1);
  free(pos);
  power=malloc(nrbins*sizeof(float));
  count=malloc(nrbins*sizeof(float));
  keffs=malloc(nrbins*sizeof(float));
  if(!power || !count || !keffs)
  {
		fprintf(stderr,"Error allocating memory for power spectrum.\n");
		exit(1);
  }
  nrbins=powerspectrum(field_dims,field,nrbins, power,count,keffs);
  free(field);
  free(power);
  free(count);
  return 0;
}
