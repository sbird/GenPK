#include "gen-pk.h"
//For CHAR_BIT
#include <limits.h>

/** \file
 * Defines a few small utility functions*/
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
           fprintf(stderr, "Usage: ./gen-pk -i filenames -j other_filenames -o outdir -c (optional) cross-corr type\n"
                           "Outputs one file per particle type, with the name PK-$TYPE-$INPUT\n"
                           "Each output file has three columns, for each bin, k_eff, P(k) and N_modes\n"
		           "if -j is specified it cross-correlates files specified under -i with files specifed under\n"
		           "-j (one output per particle type).\n"
                           "If -c is specified the code computes the cross-correlation of that particle type with \n"
		           "the CDM (type 1) (within the file specified by -i)\n");
           return;
}

std::string type_str(int type)
{
        switch(type)
        {
                case BARYON_TYPE:
                        return "by";
                case DM_TYPE:
                        return "DM";
                case NEUTRINO_TYPE:
                        return "nu";
                default:
                        return "xx";
        }
}
