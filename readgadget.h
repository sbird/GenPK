#include <stdint.h>

extern struct vector
{
  float x,y,z;
} *pos,*b;
extern int swap;

extern void swap_Nbyte(char *data,int n,int m);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
extern int64_t find_block(FILE *fd,char *label);
extern int64_t read_gadget_float(float *data,char *label,FILE *fd);
/* The final argument, if one, means it will attempt to read an old format file*/
extern int64_t read_gadget_float3(float *data,char *label,FILE *fd, int old);
extern int read_gadget_head(int *npart,double *massarr,double *time,double *redshift,double *boxsize, FILE *fd, int old);


