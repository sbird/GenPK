#include "gen-pk.h"
//For memset
#include <string.h>
#include <cmath>
//For malloc
#include <stdlib.h>

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
/** \file 
 * Defines read_fieldize() , wraps fieldize() and GadgetReader*/
int read_fieldize(float * field, GadgetReader::GSnap* snap, int type, double box, int field_dims)
{
        int64_t npart_total,toread;
        int parts=0;
        int64_t npart_stars=0;
        float *pos=NULL;
        //Initially set it to skip all types
        int skip_type=(1<<N_TYPE)-1;
        /*Stars are another type of baryons*/
        if((*snap).GetNpart(type)==0 || type==STARS_TYPE)
          return 1;
        /* Set skip_type, which should include every type *other* than the one we're interested in
         * There are N_TYPE types, so skipping all types is 2^(N_TYPES)-1 and then subtract 2^type for 
         * the one we're trying to read*/
        npart_total=(*snap).GetNpart(type);
        skip_type-=(1<<type);
        /* Zero the field*/
        for(int i=0; i<2*field_dims*field_dims*(field_dims/2+1); i++)
                field[i] = 0;
        /* Add the stars if we are using baryons. */
        if(type==BARYON_TYPE)
                npart_stars=(*snap).GetNpart(STARS_TYPE);
        npart_total+=npart_stars;
        /* Try to allocate enough memory for particle table.
         * Maximum allocation of 512**3/2 particles ~ 768MB. */
        parts=MIN(npart_total,(1<<26));
        while(!(pos=(float *)malloc(3*parts*sizeof(float)))){
                fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",3*parts*sizeof(float)/1024/1024,type);
                return 1;
        }
        toread=npart_total;
        while(toread > 0){
                if(toread < parts){
                        parts=toread;
                }
                if((*snap).GetBlock("POS ",pos,parts,0,skip_type) != parts){
                        fprintf(stderr, "Error reading particle data for type %d\n",type);
                        free(pos);
                        return 1;
                }
                /*Sanity check the data*/
                /*for(int i=0; i<npart_read; i++)
                      if(fabs(pos[3*i]) > box || fabs(pos[3*i+1]) > box || fabs(pos[3*i+2]) > box){
                              fprintf(stderr, "Part %d position is at [%e,%e,%e]! Something is wrong!\n",i,pos[3*i],pos[3*i+1],pos[3*i+2]);
                      }*/
                /* Fieldize. positions should be an array of size 3*particles
                 * (like the output of read_gadget_float3)
                 * out is an array of size [dims*dims*dims]
                 * the "extra" switch, if set to one, will assume that the output
                 * is about to be handed to an FFTW in-place routine,
                 * and set skip the last 2 places of the each row in the last dimension
                 * It is important that pos is not longer than max_int */
                fieldize(box,field_dims,field,npart_total,parts,pos, 1);
                toread-=parts;
        }
        /*Read the stars. This needs to be done in a separate GetBlocks call at the moment*/
        if(type==BARYON_TYPE){
                /*Re-allocate in case there are more stars than baryons (unlikely)*/
                if(npart_stars > parts){
                        free(pos);
                        while(!(pos=(float *)malloc(3*npart_stars*sizeof(float)))){
                                        fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",3*npart_stars*sizeof(float)/1024/1024,type);
                                        return 1;
                        }
                }
                skip_type=(1<<N_TYPE)-(1<<STARS_TYPE);
                if((*snap).GetBlock("POS ",pos,npart_stars,0,skip_type) != npart_stars){
                        fprintf(stderr, "Error reading particle data for type %d\n",type);
                        free(pos);
                        return 1;
                }
                fieldize(box,field_dims,field,npart_total,npart_stars,pos, 1);
        }
        free(pos);
        return 0;
}
