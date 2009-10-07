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
  int nfiles=1, file, segment_part;
  int64_t total_part;
  double massarr[5], time, redshift;
  double boxsize;
  float *pos, *power, *count, *keffs;
  float *field;
  FILE *fd;
  if(argc<3)
  {
			 fprintf(stderr,"Usage: NumFiles filenames\n");
			 exit(0);
  }
  //Assume single argument is a single file.
//  if(argc==2)
//         *fd=fopen(argv[1],"r");
  //Don't support old switch any more. Makes it more complicated.
//  if(argc==3 && (atoi(argc[argc-1])==1))
//      old=1;
  //Otherwise we want to read files in sequence.
  nfiles=atoi(argv[1]);
  if(nfiles < 1 || nfiles > argc-2)
  {
		 fprintf(stderr,"Filenames don't match number of files specified.\n");
		 exit(0);
  }
  /*First read the header and allocate memory.*/
  fd=fopen(argv[2],"r");
  if(!fd)
  {
		fprintf(stderr,"Error opening file %s for reading!\n", argv[2]);
		exit(1);
  }
  if(!read_gadget_head(npart, massarr, &time, &redshift,&boxsize, fd, old))
  {
		fprintf(stderr,"Error reading file header!\n");
		exit(1);
  }
  fprintf(stderr, "Boxsize=%g, npart=[%g,%g,%g,%g,%g], redshift=%g\n",boxsize,cbrt(npart[0]*nfiles),cbrt(npart[1]*nfiles),cbrt(npart[2]*nfiles),cbrt(npart[3]*nfiles),cbrt(npart[4]*nfiles),redshift);
  field_dims=(2*cbrt(npart[1]) < FIELD_DIMS ? 2*cbrt(npart[1]) : FIELD_DIMS);
  nrbins=floor(sqrt(3)*((field_dims+1.0)/2.0)+1);
  /* Here we are assuming that the particles are divided evenly amongst the files.
   * This will be the case if numfiles divides npart evenly.*/
  segment_part=npart[0]+npart[1]+npart[2]+npart[3]+npart[4];
  total_part=segment_part*nfiles;
  pos=malloc(3*(segment_part+1)*sizeof(float));
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
  fieldize(boxsize,field_dims,field,total_part, segment_part,pos, 1);
  /*Now read the rest of the files, if there are any.*/
  for(file=1; file<nfiles; file++)
  {
    fd=fopen(argv[file+2],"r");
    if(!fd)
    {
   	fprintf(stderr,"Error opening file %s for reading!\n", argv[file+2]);
   	exit(1);
    }
    if(!read_gadget_float3(pos, "POS ",fd, old))
    {
  		fprintf(stderr, "Error reading particle data\n");
  		exit(1);
    }
    //By now we should have the data.
    fclose(fd);
    fieldize(boxsize,field_dims,field,total_part,segment_part,pos, 1);
  }
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
  for(int i=0;i<nrbins;i++)
  {
    if(count[i])
      printf("%e\t%e\t%e\n",keffs[i],power[i],count[i]);
  }
  free(field);
  free(power);
  free(count);
  return 0;
}
