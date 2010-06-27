/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

#define _GNU_SOURCE
#include <string.h>
#include <limits.h>
#include "gen-pk.h"
#include <math.h>

/* In practice this means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 1024
#define PART_TYPES 6
/*Which particle type to use for stars, which are just another kind of baryon*/
#define STARS 4
int nexttwo(int);
#define MAX(x,y) ((x) > (y) ? (x) :(y))
#define MIN(x,y) ((x) < (y) ? (x) :(y))

int main(int argc, char* argv[]){
  int field_dims=0;
  int old =0, type;
  int64_t tot_npart[PART_TYPES],nrbins;
  double mass[PART_TYPES],tot_mass;
  int nfiles=1, file;
  struct gadget_header *headers;
  double boxsize, redshift;
  float *pos, *power[PART_TYPES], *count[PART_TYPES], *keffs[PART_TYPES];
/*   float *tot_power,*tot_keffs; */
  float *field;
  FILE *fd;
  if(argc<3)
  {
			 fprintf(stderr,"Usage: NumFiles filenames outdir\n");
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
  if(nfiles < 1 || nfiles > argc-3)
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
     fd=fopen(argv[2+file],"r");
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
  tot_mass=0;
  for(type=0;type<PART_TYPES;type++)
  {
          tot_npart[type]=0;
          mass[type]=headers[0].mass[type];
          tot_mass+=mass[type];
  }
  /*Assemble totals, check all the files are from the same simulation*/
  /*We can of course get totals from the header, but I don't really want to.*/
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
     {
       tot_npart[type]+=headers[file].npart[type];
       if(mass[type]!=headers[file].mass[type])
       {
         fprintf(stderr,"Error! Mass of particle %d from file 0 is %e, file %d has %e\n",type,mass[type],file,headers[file].mass[type]);
         exit(2);
       }
     }
  }
  for(type=0;type<PART_TYPES;type++)
  {
    int tmp=2*nexttwo(cbrt(tot_npart[type]));
    field_dims=MAX(field_dims, MIN(tmp, FIELD_DIMS));
  }
  nrbins=floor(sqrt(3)*((field_dims+1.0)/2.0)+1);
     fprintf(stderr, "Boxsize=%g, ",boxsize);
     fprintf(stderr, "tot_npart=[%g,%g,%g,%g,%g,%g], ",cbrt(tot_npart[0]),cbrt(tot_npart[1]),cbrt(tot_npart[2]),cbrt(tot_npart[3]),cbrt(tot_npart[4]),cbrt(tot_npart[5]));
     fprintf(stderr, "Masses=[%g %g %g %g %g %g], ",mass[0],mass[1],mass[2],mass[3],mass[4],mass[5]);
     fprintf(stderr, "redshift=%g, Ω_M=%g Ω_B=%g\n",redshift,headers[0].Omega0,mass[0]/tot_mass*headers[0].Omega0);
  /*Now read the particle data.*/
  /*Type 4 particles (stars) are treated as just another type of baryon*/
  tot_npart[0]+=tot_npart[STARS];
  tot_npart[STARS]=0;
  for(type=0; type<PART_TYPES; type++)
  {
    if(tot_npart[type]==0)
      continue;
    /* Allocating a bit more memory allows us to do in-place transforms.*/
    field=calloc(2*field_dims*field_dims*(field_dims/2+1),sizeof(float));
    if(!field)
    {
  		fprintf(stderr,"Error allocating memory for field\n");
  		exit(1);
    }
    for(file=0; file<nfiles; file++)
    {
      int npart=headers[file].npart[type];
      /*Add the stars*/
      int nstar=0;
      int offset=0;
      if(type==0)
         nstar=headers[file].npart[STARS];
      fd=fopen(argv[file+2],"r");
      if(!fd)
      {
        	fprintf(stderr,"Error opening file %s for reading!\n", argv[file+2]);
        	exit(1);
      }
      pos=malloc(3*(npart+nstar+1)*sizeof(float));
      if(!pos)
      {
    		fprintf(stderr,"Error allocating particle memory\n");
    		exit(1);
      }
      for(int i=0; i<type; i++)
              offset+=headers[file].npart[i];
      if(old) fseek(fd,sizeof(struct gadget_header)+2*sizeof(int),SEEK_CUR);
      if(read_gadget_float3(pos, "POS ",offset ,npart, fd,old) != npart)
      {
    		fprintf(stderr, "Error reading particle data\n");
    		exit(1);
      }
      /*Read stars as well. Note past this point no distinction is made between stars and baryons.*/
      if(type==0 && nstar>0){
         offset=0;
         for(int i=0; i<STARS; i++)
              offset+=headers[file].npart[i];
         if(read_gadget_float3(pos+3*npart, "POS ",offset ,nstar, fd,old) != nstar)
         {
       		fprintf(stderr, "Error reading star data\n");
    		   exit(1);
         }
      }
      //By now we should have the data.
      fclose(fd);
      /*Sanity check the data*/
      for(int i=0; i<npart+nstar; i++)
            if(fabs(pos[3*i]) > boxsize || fabs(pos[3*i+1]) > boxsize || fabs(pos[3*i+2]) > boxsize){
                    fprintf(stderr, "Part %d position is at [%e,%e,%e]! Something is wrong!\n",i,pos[3*i],pos[3*i+1],pos[3*i+2]);
            }
    /* Fieldize. positions should be an array of size 3*particles 
     * (like the output of read_gadget_float3)
     * out is an array of size [dims*dims*dims]
     * the "extra" switch, if set to one, will assume that the output 
     * is about to be handed to an FFTW in-place routine, 
     * and set skip the last 2 places of the each row in the last dimension
     */
      fieldize(boxsize,field_dims,field,tot_npart[type],npart+nstar,pos, 1);
      free(pos);
    }
    power[type]=malloc(nrbins*sizeof(float));
    count[type]=malloc(nrbins*sizeof(float));
    keffs[type]=malloc(nrbins*sizeof(float));
    if(!power[type] || !count[type] || !keffs[type])
    {
  		fprintf(stderr,"Error allocating memory for power spectrum.\n");
  		exit(1);
    }
    nrbins=powerspectrum(field_dims,field,nrbins, power[type],count[type],keffs[type]);

    free(field);
  }
  /*Print power. Note use the count from the DM particles, because 
   * they dominate the modes. I'm not sure the sample variance 
   * really decreases by a factor of two from adding a subdominant baryon component.*/

  char *filename;
  int totfil=strlen(argv[argc-1])+strlen(basename(argv[2]))+10;
  if(!(filename=malloc(totfil*sizeof(char)))){
       fprintf(stderr,"Error allocating string memory\n");
       exit(1);
  }
  /*Print out a baryon P(k) if there are any baryons*/
  if(tot_npart[0]){
     strcpy(filename,argv[argc-1]);
     strcat(filename,"/PK-by-");
     strcat(filename,basename(argv[2]));
     if(!(fd=fopen(filename, "w"))){
        fprintf(stderr,"Error opening file: %s\n",filename);
        exit(1);
     }
     for(int i=0;i<nrbins;i++)
     {
       if(count[0][i])
         fprintf(fd,"%e\t%e\t%e\n",keffs[0][i],power[0][i],count[0][i]);
     }
     fclose(fd);
  }
  /*Print out a DM P(k) if there is any DM*/
  if(tot_npart[1]){
     strcpy(filename,argv[argc-1]);
     strcat(filename,"/PK-DM-");
     strcat(filename,basename(argv[2]));
     if(!(fd=fopen(filename, "w"))){
        fprintf(stderr,"Error opening file: %s\n",filename);
        exit(1);
     }
   free(filename);
     for(int i=0;i<nrbins;i++)
     {
       if(count[1][i])
         fprintf(fd,"%e\t%e\t%e\n",keffs[1][i],power[1][i],count[1][i]);
     }
     fclose(fd);
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
  free(headers);
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
