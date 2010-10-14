#include "gen-pk.h"
//For memset
#include <string.h>
#include <cmath>

int read_fieldize(float * field, GadgetReader::GSnap* snap, int type, double box, int field_dims)
{
        int64_t npart_read;
        float *pos=NULL;
        //Initially set it to skip all types
        int skip_type=(1<<N_TYPE)-1;
        /*Stars are another type of baryons*/
        if((*snap).GetNpart(type)==0 || type==STARS_TYPE)
          return 1;
        /* Set skip_type, which should include every type *other* than the one we're interested in
         * There are N_TYPE types, so skipping all types is 2^(N_TYPES)-1 and then subtract 2^type for 
         * the one we're trying to read*/
        npart_read=(*snap).GetNpart(type);
        skip_type-=(1<<type);
        /* Add the stars if we are using baryons. */
        if(type==BARYON_TYPE){
                npart_read+=(*snap).GetNpart(STARS_TYPE);
                skip_type-=(1<<STARS_TYPE);
        }
        /*Try to allocate enough memory for particle table. If we can't, read it in chunks*/
        while(!(pos=(float *)malloc(3*(npart_read+1)*sizeof(float)))){
          	        fprintf(stderr,"Error allocating particle memory for type %d\n",type);
                  	return 1;
        }
        if((*snap).GetBlock("POS ",pos,npart_read,0,skip_type) != npart_read){
          	fprintf(stderr, "Error reading particle data for type %d\n",type);
                free(pos);
          	return 1;
        }
        /*Sanity check the data*/
        for(int i=0; i<npart_read; i++)
              if(fabs(pos[3*i]) > box || fabs(pos[3*i+1]) > box || fabs(pos[3*i+2]) > box){
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
        fieldize(box,field_dims,field,npart_read,npart_read,pos, 1);
        free(pos);
        return 0;
}
