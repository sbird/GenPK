#include "gen-pk.h"
#include <limits.h>
/* In practice this means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 1024
#define PART_TYPES 6
int nexttwo(int);

int main(char argc, char* argv[]){
  int field_dims;
  int old =0, type;
  int tot_npart[PART_TYPES],nrbins[PART_TYPES];
  int nfiles=1, file;
  struct gadget_header *headers;
  double boxsize, redshift;
  float *pos, *power[PART_TYPES], *count[PART_TYPES], *keffs[PART_TYPES];
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
  headers=malloc(nfiles*sizeof(struct gadget_header));
  if(!headers)
  {
    fprintf(stderr, "Error allocating memory for headers.\n");
    exit(1);
  }
  /*First read all the headers, allocate some memory and work out the totals.*/
  for(file=0; file<nfiles; file++)
  {
     fd=fopen(argv[2],"r");
     if(!fd)
     {
   		fprintf(stderr,"Error opening file %s for reading!\n", argv[2]);
   		exit(1);
     }
     if(!read_gadget_head(headers+file, fd, old))
     {
   		fprintf(stderr,"Error reading file header!\n");
   		exit(1);
     }
     //By now we should have the header data.
     fclose(fd);
  }
  boxsize=headers[0].BoxSize;
  redshift=headers[0].redshift;
  for(type=0;type<PART_TYPES;type++)
          tot_npart[type]=0;
  /*Assemble totals, check all the files are from the same simulation*/
  for(file=0;file<nfiles;file++)
  {
     if(boxsize!=headers[file].BoxSize)
     {
       fprintf(stderr,"Error! Box size from file 0 is %e, while file %d has %e\n",boxsize,file,headers[file].BoxSize);
       exit(2);
     }
     if(redshift!=headers[file].redshift)
     {
       fprintf(stderr,"Error! Redshift from file 0 is %e, while file %d has %e\n",redshift,file,headers[file].redshift);
       exit(2);
     }
     for(type=0;type<PART_TYPES;type++)
       tot_npart[type]+=headers[file].npart[type];
  }
     fprintf(stderr, "Boxsize=%g, tot_npart=[%g,%g,%g,%g,%g,%g], redshift=%g\n",boxsize,cbrt(tot_npart[0]),cbrt(tot_npart[1]),cbrt(tot_npart[2]),cbrt(tot_npart[3]),cbrt(tot_npart[4]),cbrt(tot_npart[5]),redshift);
  /*Now read the particle data.*/
  for(type=0; type<PART_TYPES; type++)
  {
    if(tot_npart[type]==0)
      continue;
    /* Definitely a faster bitshifting way to do this. */
    int tmp=2*nexttwo(cbrt(tot_npart[type]));
    field_dims=(tmp < FIELD_DIMS ? tmp : FIELD_DIMS);
    nrbins[type]=floor(sqrt(3)*((field_dims+1.0)/2.0)+1);
    /* Allocating a bit more memory allows us to do in-place transforms.*/
    field=malloc(2*field_dims*field_dims*(field_dims/2+1)*sizeof(float));
    if(!field)
    {
  		fprintf(stderr,"Error allocating memory for field\n");
  		exit(1);
    }
    for(file=0; file<nfiles; file++)
    {
      int npart=headers[file].npart[type];
      fd=fopen(argv[file+2],"r");
      if(!fd)
      {
        	fprintf(stderr,"Error opening file %s for reading!\n", argv[file+2]);
        	exit(1);
      }
      pos=malloc(3*(npart+1)*sizeof(float));
      if(!pos)
      {
    		fprintf(stderr,"Error allocating particle memory\n");
    		exit(1);
      }
      if(read_gadget_float3(pos, "POS ",(type<0 ? 0 : headers[file].npart[type-1]) ,npart, fd,old) != npart)
      {
    		fprintf(stderr, "Error reading particle data\n");
    		exit(1);
      }
      //By now we should have the data.
      fclose(fd);
    /* Fieldize. positions should be an array of size 3*particles 
     * (like the output of read_gadget_float3)
     * out is an array of size [dims*dims*dims]
     * the "extra" switch, if set to one, will assume that the output 
     * is about to be handed to an FFTW in-place routine, 
     * and set skip the last 2 places of the each row in the last dimension
     */
      fieldize(boxsize,field_dims,field,tot_npart[type],npart,pos, 1);
      free(pos);
    }
    power[type]=malloc(nrbins[type]*sizeof(float));
    count[type]=malloc(nrbins[type]*sizeof(float));
    keffs[type]=malloc(nrbins[type]*sizeof(float));
    if(!power[type] || !count[type] || !keffs[type])
    {
  		fprintf(stderr,"Error allocating memory for power spectrum.\n");
  		exit(1);
    }
    nrbins[type]=powerspectrum(field_dims,field,nrbins[type], power[type],count[type],keffs[type]);
    free(field);
    fprintf(stderr, "Type %d done\n",type);
  }
  /*Only print out DM particles for now*/
  for(int i=0;i<nrbins[0];i++)
  {
    if(count[0][i])
      printf("%e\t%e\t%e\n",keffs[0][i],power[0][i],count[0][i]);
  }
  for(type=0; type<PART_TYPES; type++)
  {
    if(tot_npart[type])
    {
     free(power[type]);
     free(count[type]);
     free(keffs[type]);
    }
  }
  return 0;
}

/*Returns the maximum value of an array of size size*/
/*int maxarr(int *arr, int size)
{
   int max=*arr;
   while(arr<arr+size)
   {
      max=(max > *(++arr) ? max : *arr);
   }
   return max;
}*/

/*Returns the next power of two. Stolen from wikipedia.*/
int nexttwo(int n)
{
    int i;
    n--;
    for(i=1;i<sizeof(int)*CHAR_BIT; i<<=1)
       n |= n>>i;
    return ++n; 
}
