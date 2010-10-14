/* Copyright (c) 2009,2010 Simeon Bird <spb41@cam.ac.uk>
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

#include "gadgetreader.hpp"
#include "gen-pk.h"
#include <math.h>
//For getopt
#include <unistd.h>
//For CHAR_BIT
#include <limits.h>
//For omp_get_num_procs
#include <omp.h>
#include <string>
//For memset
#include <string.h>
#include <stdlib.h>

/* In practice this means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 1024
int nexttwo(int);
int print_pk(std::string filename, int nrbins, float * keffs, float * power, int* count);
// #define MAX(x,y) ((x) > (y) ? (x) :(y))
// #define MIN(x,y) ((x) < (y) ? (x) :(y))

using namespace GadgetReader;
using namespace std;

int main(int argc, char* argv[]){
  int nrbins,field_dims=0,type;
  float *pos,*field, *power[N_TYPE], *keffs[N_TYPE];
  int *count[N_TYPE]; 
  string infiles(""),outdir("");
  char c;
  gadget_header head;
  fftwf_plan pl;
  fftwf_complex *outfield;
  while((c = getopt(argc, argv, "i:o:h")) !=-1){
    switch(c){
        case 'o':
           outdir=static_cast<string>(optarg);
           break;
        case 'i':
           infiles=static_cast<string>(optarg);
           break;
        case 'h':
        default:
           help();
           return 0;
      }
  }
  //Open the snapshot
  GSnap snap(infiles);
  if(outdir.empty() || snap.GetNumFiles() < 1){
          help();
          return 0;
  }
  //Work out how large a field we need
  for(type=0;type<N_TYPE;type++){
    int tmp=2*nexttwo(cbrt(snap.GetNpart(type)));
    field_dims=std::max(field_dims, std::min(tmp, FIELD_DIMS));
  }
  //Get the header and print out some useful things
  head=snap.GetHeader();
  nrbins=floor(sqrt(3)*((field_dims+1.0)/2.0)+1);
  fprintf(stderr, "Boxsize=%g, ",head.BoxSize);
  fprintf(stderr, "redshift=%g, Ω_M=%g\n",head.redshift,head.Omega0);
  fprintf(stderr, "NPart=(%g,%g,%g,%g,%g,%g)**3\n",cbrt(snap.GetNpart(0)),cbrt(snap.GetNpart(1)),cbrt(snap.GetNpart(2)),cbrt(snap.GetNpart(3)),cbrt(snap.GetNpart(4)),cbrt(snap.GetNpart(5)));
  fprintf(stderr, "Masses=[%g %g %g %g %g %g]\n",head.mass[0],head.mass[1],head.mass[2],head.mass[3],head.mass[4],head.mass[5]);
  /*Now make a power spectrum for each particle type*/
  /* Allocating a bit more memory allows us to do in-place transforms.*/
  if(!(field=(float *)fftwf_malloc(2*field_dims*field_dims*(field_dims/2+1)*sizeof(float)))){
  	fprintf(stderr,"Error allocating memory for field\n");
  	exit(1);
  }
  outfield=(fftwf_complex *) &field[0];
  if(!fftwf_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  return 0;
  }
  fftwf_plan_with_nthreads(omp_get_num_procs());
  pl=fftwf_plan_dft_r2c_3d(field_dims,field_dims,field_dims,&field[0],outfield, FFTW_ESTIMATE);
  for(type=0; type<N_TYPE; type++)
  {
        int64_t npart_read;
        //Initially set it to skip all types
        int skip_type=(1<<N_TYPE)-1;
        /*Stars are another type of baryons*/
        if(snap.GetNpart(type)==0 || type==STARS_TYPE)
          continue;
        /* Set skip_type, which should include every type *other* than the one we're interested in
         * There are N_TYPE types, so skipping all types is 2^(N_TYPES)-1 and then subtract 2^type for 
         * the one we're trying to read*/
        npart_read=snap.GetNpart(type);
        skip_type-=(1<<type);
        /* Add the stars if we are using baryons. */
        if(type==BARYON_TYPE){
                npart_read+=snap.GetNpart(STARS_TYPE);
                skip_type-=(1<<STARS_TYPE);
        }
        /*Try to allocate enough memory for particle table. If we can't, read it in chunks*/
        while(!(pos=(float *)malloc(3*(npart_read+1)*sizeof(float)))){
          	        fprintf(stderr,"Error allocating particle memory for type %d\n",type);
                  	continue;
        }
        if(snap.GetBlock("POS ",pos,npart_read,0,skip_type) != npart_read){
          	fprintf(stderr, "Error reading particle data for type %d\n",type);
                free(pos);
          	continue;
        }
        /*Sanity check the data*/
        for(int i=0; i<npart_read; i++)
              if(fabs(pos[3*i]) > head.BoxSize || fabs(pos[3*i+1]) > head.BoxSize || fabs(pos[3*i+2]) > head.BoxSize){
                      fprintf(stderr, "Part %d position is at [%e,%e,%e]! Something is wrong!\n",i,pos[3*i],pos[3*i+1],pos[3*i+2]);
              }
        /* Fieldize. positions should be an array of size 3*particles 
         * (like the output of read_gadget_float3)
         * out is an array of size [dims*dims*dims]
         * the "extra" switch, if set to one, will assume that the output 
         * is about to be handed to an FFTW in-place routine, 
         * and set skip the last 2 places of the each row in the last dimension
         */
        memset(field,0,2*field_dims*field_dims*(field_dims/2+1));
        fieldize(head.BoxSize,field_dims,field,snap.GetNpart(type),npart_read,pos, 1);
        free(pos);
        power[type]=(float *) malloc(nrbins*sizeof(float));
        count[type]=(int *) malloc(nrbins*sizeof(int));
        keffs[type]=(float *) malloc(nrbins*sizeof(float));
        if(!power[type] || !count[type] || !keffs[type]){
      		fprintf(stderr,"Error allocating memory for power spectrum.\n");
                continue;
        }
        nrbins=powerspectrum(field_dims,&pl,outfield,nrbins, power[type],count[type],keffs[type]);
  }
  fftwf_free(field);
  fftwf_destroy_plan(pl);
  string filename=outdir;
  size_t last=infiles.find_last_of("/\\");
  /*Print power. Note use the count from the DM particles, because 
   * they dominate the modes. I'm not sure the sample variance 
   * really decreases by a factor of two from adding a subdominant baryon component.*/
  /*Print out a baryon P(k) if there are any baryons*/
  if(snap.GetNpart(BARYON_TYPE)){
     filename+="/PK-by-"+infiles.substr(last+1);
     print_pk(filename,nrbins,keffs[BARYON_TYPE],power[BARYON_TYPE],count[DM_TYPE]);
  }
  /*Print out a DM P(k) if there is any DM*/
  if(snap.GetNpart(DM_TYPE)){
     filename=outdir;
     filename+="/PK-DM-"+infiles.substr(last+1);
     print_pk(filename,nrbins,keffs[DM_TYPE],power[DM_TYPE],count[DM_TYPE]);
  }
  for(type=0; type<N_TYPE; type++){
    if(type != STARS_TYPE && snap.GetNpart(type)){
         free(power[type]);
         free(count[type]);
         free(keffs[type]);
    }
  }
  return 0;
}

int print_pk(std::string filename, int nrbins, float * keffs, float * power, int * count)
{
  FILE *fd;
  if(!(fd=fopen(filename.c_str(), "w"))){
     fprintf(stderr,"Error opening file: %s\n",filename.c_str());
     return 0;
  }
  for(int i=0;i<nrbins;i++)
  {
    if(count[i])
      fprintf(fd,"%e\t%e\t%d\n",keffs[i],power[i],count[i]);
  }
  fclose(fd);
  return nrbins;
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
    unsigned int i;
    n--;
    for(i=1;i<sizeof(int)*CHAR_BIT; i<<=1)
       n |= n>>i;
    return ++n; 
}

void help()
{
           fprintf(stderr, "Usage: ./gen-pk -i filenames -o outdir\n");
           return;
}
