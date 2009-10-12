#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>

extern struct vector
{
  float x,y,z;
} *pos,*b;
extern int swap;

struct gadget_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      numfiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes
 */
};

extern void swap_Nbyte(char *data,int n,int m);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
extern int64_t find_block(FILE *fd,char *label);
extern int64_t read_gadget_float(float *data,char *label,FILE *fd);
/* The final argument, if one, means it will attempt to read an old format file*/
extern int64_t read_gadget_float3(float *data,char *label,int offset, int read, FILE *fd, int old);
extern int read_gadget_head(struct gadget_header *out_header, FILE *fd, int old);


