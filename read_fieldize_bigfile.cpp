#include "gen-pk.h"
#ifndef NOBIGFILE
#include "bigfile/src/bigfile.h"

#ifdef DOUBLE_PRECISION
    #define H5_DATASET_READ H5LTread_dataset_double
#else
    #define H5_DATASET_READ H5LTread_dataset_float
#endif
/*Returns true if an be opened as a BigFile*/
int is_bigfile(const char * infile)
{
  BigFile bf = {0};
  if(0 != big_file_open(&bf, infile))
      return 0;
  big_file_close(&bf);
  return 1;
}

/* this routine loads header data from the first file of an HDF5 snapshot.*/
int load_bigfile_header(const char *fname, double  *atime, double *redshift, double *box100, double *h100, int64_t *npart_all, double * mass, double *Omega0)
{
  double OmegaLambda;
  BigFile bf = {0};
  if(0 != big_file_open(&bf, fname)) {
      fprintf(stderr,"Failed to open snapshot at %s:%s\n", fname,
                  big_file_get_error_message());
  }
  BigBlock bh = {0};
  if(0 != big_file_open_block(&bf, &bh, "header")) {
      fprintf(stderr,"Failed to create block at %s:%s\n", "header",
                  big_file_get_error_message());
  }
  if(
  (0 != big_block_get_attr(&bh, "TotNumPart", npart_all, "u8", N_TYPE)) ||
  (0 != big_block_get_attr(&bh, "MassTable", mass, "f8", N_TYPE)) ||
  (0 != big_block_get_attr(&bh, "Time", atime, "f8", 1)) ||
  (0 != big_block_get_attr(&bh, "Redshift", redshift, "f8", 1)) ||
  (0 != big_block_get_attr(&bh, "HubbleParam", h100, "f8", 1)) ||
  (0 != big_block_get_attr(&bh, "Omega0", Omega0, "f8", 1)) ||
  (0 != big_block_get_attr(&bh, "OmegaL", &OmegaLambda, "f8", 1)) ||
  (0 != big_block_get_attr(&bh, "BoxSize", box100, "f8", 1))) {
      fprintf(stderr,"Failed to read attr: %s\n",
                  big_file_get_error_message());
  }
  if(big_block_close(&bh) ||
        big_file_close(&bf)) {
      fprintf(stderr,"Failed to close block or file: %s\n",
                  big_file_get_error_message());
  }
  printf("NumPart=[%ld,%ld,%ld,%ld,%ld,%ld], ",npart_all[0],npart_all[1],npart_all[2],npart_all[3],npart_all[4],npart_all[5]);
  printf("Masses=[%g %g %g %g %g %g], ",mass[0],mass[1],mass[2],mass[3],mass[4],mass[5]);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",(*redshift),*Omega0,OmegaLambda);
  printf("Expansion factor = %f\n",(*atime));
  printf("Hubble = %g Box=%g \n",(*h100),(*box100));
  return 0;
}
  
/* This routine loads particle data from a bigfile snapshot set. */
int read_fieldize_bigfile(GENFLOAT * field, const char *fname, int type, double box, int field_dims, double *total_mass, int64_t* npart_all, double * mass, double Omega0)
{
  char name[32];
  BigFile bf = {0};
  if(0 != big_file_open(&bf, fname)) {
      fprintf(stderr,"Failed to open snapshot at %s:%s\n", fname,
                  big_file_get_error_message());
      return 1;
  }
  BigBlock bb = {0};
  snprintf(name,32,"%d/Position",type);
  if(0 != big_file_open_block(&bf, &bb, name)) {
      fprintf(stderr,"Failed to create block at %s:%s\n", "Position",
                  big_file_get_error_message());
      return 1;
  }
  BigArray pos = {0};
  /*Open particle data*/
  if(0 != big_block_read_simple(&bb, 0, npart_all[type], &pos,"f8")) {
      fprintf(stderr,"Failed to read from block: %s\n", big_file_get_error_message());
      return 1;
  }
  big_block_close(&bb);

  float * positions = (float *) malloc(3*npart_all[type]*sizeof(float));
  for(int i=0; i<3*npart_all[type]; i++)
      positions[i] = ((double *)pos.data)[i];
  free(pos.data);
  BigArray massarray = {0};
  /* Load particle masses, if present  */
   if(mass[type] == 0){
        double total_mass_this_file=0;
        snprintf(name,32,"%d/Masses",type);
        if(0 != big_block_read_simple(&bb, 0, npart_all[type], &massarray,"f4")) {
            fprintf(stderr,"Failed to read from block: %s\n", big_file_get_error_message());
            return 1;
        }
        for(int i = 0; i<npart_all[type]; i++)
            total_mass_this_file += ((float *)massarray.data)[i];
        *total_mass += total_mass_this_file;
  }
  if(big_block_close(&bb) ||
        big_file_close(&bf)) {
      fprintf(stderr,"Failed to close block or file: %s\n",
                  big_file_get_error_message());
  }
  //Do the final summation here to avoid fp roundoff
  *total_mass += mass[type]*npart_all[type];
  fieldize(box,field_dims,field,npart_all[type],positions, (float *)massarray.data, mass[type], 1);
  free(positions);
  free(massarray.data);
  return 0;
}
#endif //NOBIGFILE

