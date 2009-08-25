#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readgadget.h"
#define int4bytes int
/*--------- comment/uncomment to remove/enable DEBUG outputs ------------------*/
/*
#define MY_DEBUG
*/

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Low Level Routines -----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

int4bytes blksize,swap=0;
#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte((char*)&blksize,1,4);}

/*---------------------- Basic routine to read data from a file ---------------*/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) !\n");
      exit(3);
    }
  return nread;
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data,int n,int m)
{
  int i,j;
  char old_data[16];

  if(swap>0)
    {
      for(j=0;j<n;j++)
	{
          memcpy(&old_data[0],&data[j*m],m);
          for(i=0;i<m;i++)
            {
              data[j*m+i]=old_data[m-i-1];
	    }
	}
    }
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
int find_block(FILE *fd,char *label)
{
  int4bytes blocksize=0;
  char blocklabel[5]={"    "};

  rewind(fd);

  while(!feof(fd) && blocksize == 0)
    {
       SKIP;
       if(blksize == 134217728)
	 {
#ifdef MY_DEBUG
           printf("Enable ENDIAN swapping !\n");
#endif
           swap=1-swap;
           swap_Nbyte((char*)&blksize,1,4);
	 }
       if(blksize != 8)
         {
	   printf("incorrect format (blksize=%d)!\n",blksize);
           exit(1);
         }
       else
         {
           my_fread(blocklabel, 4*sizeof(char), 1, fd);
           my_fread(&blocksize, sizeof(int4bytes), 1, fd);
           swap_Nbyte((char*)&blocksize,1,4);
#ifdef MY_DEBUG
           printf("Found Block <%s> with %d bytes\n",blocklabel,blocksize);
#endif
           SKIP;
	   if(strcmp(label,blocklabel)!=0)
	     { 
                fseek(fd,blocksize,1);
                blocksize=0;
	     }
         }
    }
  return(blocksize-8);
}


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- High Level Routines ----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read the header information ---------------*/
/*-------- int *npart:    List of Particle numbers for spezies [0..5] ---------*/
/*-------- int *massarr:  List of masses for spezies [0..5] -------------------*/
/*-------- int *time:     Time of snapshot ------------------------------------*/
/*-------- int *redshift:  Redshift of snapshot --------------------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns number of read bytes ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_head(int *npart,double *massarr,double *time,double *redshift,FILE *fd)
{
  int blocksize,dummysize,i;

  blocksize = find_block(fd,"HEAD");
  if(blocksize <= 0)
    {
      printf("Block <%s> not fond !\n","HEAD");
      exit(5);
    }
  else
    {
       dummysize=blocksize - 6 * sizeof(int) - 8 * sizeof(double);
       SKIP;
       my_fread(npart,6*sizeof(int), 1, fd);        swap_Nbyte((char*)npart,6,4);
       my_fread(massarr,6*sizeof(double), 1, fd);   swap_Nbyte((char*)massarr,6,8);
       my_fread(time,sizeof(double), 1, fd);        swap_Nbyte((char*)time,1,8);
       my_fread(redshift,sizeof(double), 1, fd);    swap_Nbyte((char*)redshift,1,8);
		 //flag_sfr,flag_feedback,partTotal, $
         //                flag_cooling,num_files,BoxSize,Omega0,OmegaLambda,HubbleParam
       fseek(fd,dummysize,1);
       SKIP;
    }
  return(blocksize);
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 1D float array ---------------------*/
/*-------- float *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_float(float *data,char *label,FILE *fd)
{
  int blocksize;

  blocksize = find_block(fd,label);
  if(blocksize <= 0)
    {
      printf("Block <%s> not fond !\n",label);
      exit(5);
    }
  else
    {
#ifdef MY_DEBUG
       printf("Reading %d bytes of data from <%s>...\n",blocksize,label);
#endif
       SKIP;
       my_fread(data,blocksize, 1, fd);
       swap_Nbyte((char*)data,blocksize/sizeof(float),4);
       SKIP;
    }
  return(blocksize/sizeof(float));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 3D float array ---------------------*/
/*-------- float *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
int read_gadget_float3(float *data,char *label,FILE *fd)
{
  int blocksize,i;

  blocksize = find_block(fd,label);
  if(blocksize <= 0)
    {
      printf("Block <%s> not fond !\n",label);
      exit(5);
    }
  else
    {
#ifdef MY_DEBUG
       printf("Reding %d bytes of data from <%s>...\n",blocksize,label);
#endif
       SKIP;
       my_fread(data,blocksize, 1, fd);
       swap_Nbyte((char*)data,blocksize/sizeof(float),4);
       SKIP;
    }
  return(blocksize/sizeof(float)/3);
}


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- Small Example HowToUse -------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

#ifdef MAKE_READGADGET_TEST
int main(int argc, char **argv)
{
  FILE *fd = 0;
  char filename[265];
  int i,n,ntot;
  int npart[6];
  double masstab[6],redshift,time;
  float *rho;
  struct vector *pos,*b;

  if(argc >= 2)
    {
      strcpy(filename,argv[1]);
      if(!(fd = fopen(filename,"r")))
        {
           printf("Cant open file <%s> !\n",filename);
           exit(2);
        }  
      else
        {
	  /*----------- READ HEADER TO GET GLOBAL PROPERTIES -------------*/
 	   n = read_gadget_head(npart,masstab,&time,&redshift,fd);

           ntot=0;
           for(i=0;i<6;i++)
	     {
	       printf("PartSpezies %d, anz=%d, masstab=%f\n",i,npart[i],masstab[i]);
               ntot += npart[i];
	     }
           printf("Time of snapshot=%f, z=%f, ntot=%d\n",time,redshift,ntot);

	   /*---------- ALLOCATE MEMORY ---------------------------------*/
           rho=malloc(npart[0]*sizeof(float));
           b=malloc(3*npart[0]*sizeof(float));
           pos=malloc(3*ntot*sizeof(float));

	   /*---------- READ DATA BLOCKS --------------------------------*/
	   n = read_gadget_float(rho,"RHO ",fd);
	   n = read_gadget_float3((float*)pos,"POS ",fd);
	   n = read_gadget_float3((float*)b,"BFLD",fd);
	   /*
           for(i=0;i<npart[0];i++)
	     {
	        printf("%d: (%f %f %f) =  %f\n",i,pos[i].x,pos[i].y,pos[i].z,rho[i]);
	     }
	   */
	   /*---------- CLOSE FILE AND FREE DATA ------------------------*/
	   fclose(fd);

           free(rho);
           free(b);
           free(pos);
	}
    }
  else
    {
      printf("Please give a filename ...\n");
      exit(4);
    }
  exit(0);
} 
#endif










