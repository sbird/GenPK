#include "gen-pk.h"
//For memset
#include <string.h>
#include <cmath>
//For malloc
#include <stdlib.h>
#ifndef NOHDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#include <sstream>
#endif

#define MIN(a,b) ( (a) < (b) ? (a) : (b))
/** \file
 * Defines read_fieldize() , wraps fieldize() and GadgetReader */

/** read_fieldize: reads particles from a Gadget format particle catalogue and places them on a grid, using cloud-in-cell*/
int read_fieldize(GENFLOAT * field, GadgetReader::GSnap* snap, int type, double box, int field_dims, double * total_mass)
{
        int64_t npart_total,toread;
        double mass;
        int parts=0;
        int64_t read=0;
        GENPK_FLOAT_TYPE *pos=NULL;
        GENPK_FLOAT_TYPE *masses=NULL;
        //Initially set it to skip all types
        int skip_type=(1<<N_TYPE)-1;
        /*Stars are another type of baryons*/
        if((*snap).GetNpart(type)==0)
          return 1;
        /* Set skip_type, which should include every type *other* than the one we're interested in
         * There are N_TYPE types, so skipping all types is 2^(N_TYPES)-1 and then subtract 2^type for
         * the one we're trying to read*/
        npart_total=(*snap).GetNpart(type);
        gadget_header head = (*snap).GetHeader();
        mass = head.mass[type];
        skip_type-=(1<<type);
        int64_t partlen = snap->GetBlockSize("POS ",-1) / snap->GetBlockParts("POS ");
        if ( partlen != 3*sizeof(GENPK_FLOAT_TYPE) ) {
            fprintf(stderr, "The pos array uses %ld bytes per particle, instead of %lu.\n You probably want to recompile a different DOUBLE_PRECISION_SNAP define.\n", partlen, 3*sizeof(GENPK_FLOAT_TYPE));
            exit(1);
        }
        /* Try to allocate enough memory for particle table.
         * Maximum allocation of 512**3/2 particles ~ 768MB. */
        parts=MIN(npart_total,(1<<26));
        if(!(pos=(GENPK_FLOAT_TYPE *)malloc(3*parts*sizeof(GENPK_FLOAT_TYPE)))){
                fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",3*parts*sizeof(GENPK_FLOAT_TYPE)/1024/1024,type);
                return 1;
        }
        toread=npart_total;
        while(toread > 0){
                if(toread < parts){
                        parts=toread;
                }
                if((*snap).GetBlock("POS ",pos,parts,read,skip_type) != parts){
                        fprintf(stderr, "Error reading particle data for type %d\n",type);
                        free(pos);
                        return 1;
                }
                /* Load particle masses, if present  */
                if(mass == 0){
                        double total_mass_this_file=0;
                        if(!(masses=(GENPK_FLOAT_TYPE *)malloc(parts*sizeof(GENPK_FLOAT_TYPE)))){
                            fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",parts*sizeof(GENPK_FLOAT_TYPE)/1024/1024,type);
                            return 1;
                        }
                        if((*snap).GetBlock("MASS",masses,parts,read,skip_type) != parts){
                            fprintf(stderr, "Error reading mass data for type %d\n",type);
                            free(masses);
                            return 1;
                        }
                        for(int i = 0; i<parts; i++)
                            total_mass_this_file += masses[i];
                        *total_mass += total_mass_this_file;
                }
                //Do the final summation here to avoid fp roundoff
                *total_mass += mass*parts;
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
                fieldize(box,field_dims,field,parts,pos, masses, mass, 1);
                toread-=parts;
                read+=parts;
        }
        *total_mass += 1;
        free(pos);
        return 0;
}

#ifndef NOHDF5

std::string find_first_hdf_file(const std::string& infname)
{
  /*Switch off error handling so that we can check whether a
   * file is HDF5 */
  /* Save old error handler */
  hid_t error_stack=0;
  herr_t (*old_func)(hid_t, void*);
  void *old_client_data;
  H5Eget_auto(error_stack, &old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(error_stack, NULL, NULL);
  std::string fname = infname;

  /*Were we handed an HDF5 file?*/
  if(H5Fis_hdf5(fname.c_str()) <= 0){
     /*If we weren't, were we handed an HDF5 file without the suffix?*/
     fname = infname+std::string(".0.hdf5");
     if (H5Fis_hdf5(fname.c_str()) <= 0) {
        fname = infname+std::string(".hdf5");
        if (H5Fis_hdf5(fname.c_str()) <= 0) {
            fname = std::string();
        }
     }
  }

  /* Restore previous error handler */
  H5Eset_auto(error_stack, old_func, old_client_data);
  return fname;
}

/*Open a file for reading to check it exists*/
int file_readable(const char * filename)
{
     FILE * file;
     if ((file = fopen(filename, "r"))){
          fclose(file);
          return 1;
     }
     return 0;
}

std::vector<std::string> find_hdf_set(const std::string& infname)
{
    unsigned i_fileno=0;
    std::vector<std::string> files;
    std::string fname = find_first_hdf_file(infname);
    if ( fname.empty() )
        return files;
    //Add first file
    files.push_back(fname);
    //Find position of file index in order to replace it
	i_fileno = fname.find(".0.hdf5")+1;
    //If can't find the file index, not a set but a single file.
    if(i_fileno == std::string::npos)
        return files;
    //Find filename, starting at 1
    for(int fileno = 1; fileno < 1000; fileno++) {
        std::string ffname;
        std::ostringstream convert;
        convert<<fileno;
        ffname = fname.replace(i_fileno, 1 + (fileno > 10) + (fileno > 100), convert.str());
        /*If we ran out of files, we're done*/
        if(!(file_readable(ffname.c_str()) && H5Fis_hdf5(ffname.c_str()) > 0))
                break;
        files.push_back(ffname);
    }
    return files;
}

#ifdef DOUBLE_PRECISION_SNAP
    #define H5_DATASET_READ H5LTread_dataset_double
#else
    #define H5_DATASET_READ H5LTread_dataset_float
#endif

/*Routine that is a wrapper around HDF5's dataset access routines to do error checking. Returns the length on success, 0 on failure.*/
hsize_t get_single_dataset(const char *name, GENPK_FLOAT_TYPE * data_ptr,  hsize_t data_length, hid_t * hdf_group,int fileno)
{
          int rank;
          hsize_t vlength;
          size_t type_size;
          H5T_class_t class_id;
          if (H5LTget_dataset_ndims(*hdf_group, name, &rank) < 0 || rank != 1){
             fprintf(stderr, "File %d: Rank of %s is %d !=1\n",fileno,name, rank);
             return 0;
          }
          H5LTget_dataset_info(*hdf_group, name, &vlength, &class_id, &type_size);
          if(type_size != 4 || class_id != H5T_FLOAT  || vlength > data_length || H5_DATASET_READ(*hdf_group, name, data_ptr) < 0 ){
              fprintf(stderr, "File %d: Failed reading %s (%lu)\n",fileno,name, (uint64_t)vlength);
              return 0;
          }
          return vlength;
}

/*A similar wrapper around HDF5's dataset access routines to do error checking. Returns the length on success, 0 on failure.*/
hsize_t get_triple_dataset(const char *name, GENPK_FLOAT_TYPE * data_ptr, hsize_t data_length, hid_t * hdf_group,int fileno)
{
          int rank;
          hsize_t vlength[2];
          size_t type_size;
          H5T_class_t class_id;
          if (H5LTget_dataset_ndims(*hdf_group, name, &rank) < 0 || rank != 2){
             fprintf(stderr, "File %d: Rank of %s is %d !=2\n",fileno,name, rank);
             return 0;
          }
          H5LTget_dataset_info(*hdf_group, name, &vlength[0], &class_id, &type_size);
          if(type_size != 4 || class_id != H5T_FLOAT || vlength[1] != 3 || vlength[0] > data_length || H5_DATASET_READ(*hdf_group, name, data_ptr) < 0 ){
              fprintf(stderr, "File %d: Failed reading %s (%lu)\n",fileno,name, (uint64_t)vlength[0]);
              return 0;
          }
          return vlength[0];
}

/* this routine loads header data from the first file of an HDF5 snapshot.*/
int load_hdf5_header(const char *ffname, double  *atime, double *redshift, double *box100, double *h100, int64_t *npart_all, double * mass)
{
  int i;
  unsigned int npart[N_TYPE];
  int flag_cooling;
  double Omega0, OmegaLambda;
  hid_t hdf_group,hdf_file;
  hdf_file=H5Fopen(ffname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(hdf_file < 0){
        return -1;
  }
  if ( (hdf_group=H5Gopen(hdf_file,"/Header",H5P_DEFAULT)) < 0) {
        H5Fclose(hdf_file);
        return -1;
  }
  /* Read some header functions */

  if(H5LTget_attribute_double(hdf_group,".","Time",atime) ||
     H5LTget_attribute_double(hdf_group,".","Redshift", redshift) ||
     H5LTget_attribute_double(hdf_group,".","BoxSize", box100) ||
     H5LTget_attribute_double(hdf_group,".","HubbleParam", h100) ||
     H5LTget_attribute_double(hdf_group,".","Omega0", &Omega0) ||
     H5LTget_attribute_double(hdf_group,".","OmegaLambda", &OmegaLambda) ||
     H5LTget_attribute_int(hdf_group,".","Flag_Cooling",&flag_cooling)){
          fprintf(stderr,"Failed to read some header value\n");
      H5Gclose(hdf_group);
      H5Fclose(hdf_file);
      return -1;
  }
  /*Get the total number of particles*/
  H5LTget_attribute(hdf_group,".","NumPart_Total",H5T_NATIVE_INT, &npart);
  for(i = 0; i< N_TYPE; i++)
          npart_all[i]=npart[i];
  H5LTget_attribute(hdf_group,".","NumPart_Total_HighWord",H5T_NATIVE_INT, &npart);
  for(i = 0; i< N_TYPE; i++)
          npart_all[i]+=(1L<<32)*npart[i];
  H5LTget_attribute(hdf_group,".","MassTable",H5T_NATIVE_DOUBLE, mass);

  /*Close header*/
  H5Gclose(hdf_group);
  H5Fclose(hdf_file);

  printf("NumPart=[%ld,%ld,%ld,%ld,%ld,%ld], ",npart_all[0],npart_all[1],npart_all[2],npart_all[3],npart_all[4],npart_all[5]);
  printf("Masses=[%g %g %g %g %g %g], ",mass[0],mass[1],mass[2],mass[3],mass[4],mass[5]);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",(*redshift),Omega0,OmegaLambda);
  printf("Expansion factor = %f\n",(*atime));
  printf("Hubble = %g Box=%g \n",(*h100),(*box100));
  return 0;
}

/* This routine loads particle data from a single HDF5 snapshot file.
 * A snapshot may be distributed into multiple files. */
int read_fieldize_hdf5(GENFLOAT * field, const char *ffname, int type, double box, int field_dims, double * total_mass, int fileno)
{
  int npart[N_TYPE];
  double mass[N_TYPE];
  GENPK_FLOAT_TYPE *pos=NULL;
  GENPK_FLOAT_TYPE *masses=NULL;
  int ret = 0;
  char name[16];
  double Omega0;
  hid_t hdf_group,hdf_file;
  hsize_t length;
  hdf_file=H5Fopen(ffname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(hdf_file < 0){
        return -1;
  }
  if ( (hdf_group=H5Gopen(hdf_file,"/Header",H5P_DEFAULT)) < 0) {
        H5Fclose(hdf_file);
        return -1;
  }
  if( H5LTget_attribute(hdf_group,".","NumPart_ThisFile",H5T_NATIVE_INT, &npart) ||
      H5LTget_attribute_double(hdf_group,".","Omega0", &Omega0) ||
      H5LTget_attribute(hdf_group,".","MassTable",H5T_NATIVE_DOUBLE, mass) ) {
      fprintf(stderr,"Failed to read some header value\n");
      H5Gclose(hdf_group);
      H5Fclose(hdf_file);
      return -1;
  }
  if(!(pos=(GENPK_FLOAT_TYPE *)malloc(3*npart[type]*sizeof(GENPK_FLOAT_TYPE)))){
          fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",3*npart[type]*sizeof(GENPK_FLOAT_TYPE)/1024/1024,type);
          return -1;
  }
  H5Gclose(hdf_group);
  /*If there are none of this particle in this file, do nothing.*/
  if(npart[type] == 0)
      return 0;
  /*Open particle data*/
  snprintf(name,16,"/PartType%d",type);

  if ( (hdf_group=H5Gopen(hdf_file,name,H5P_DEFAULT)) < 0) {
        H5Fclose(hdf_file);
        return -1;
  }

  /* Read position and velocity*/
  length = get_triple_dataset("Coordinates",pos,npart[type],&hdf_group,fileno);
  if(length == 0) {
          ret = -1;
          goto exit;
  }
  /* Load particle masses, if present  */
   if(mass[type] == 0){
       double total_mass_this_file=0;
         if(!(masses=(GENPK_FLOAT_TYPE *)malloc(npart[type]*sizeof(GENPK_FLOAT_TYPE)))){
          fprintf(stderr,"Error allocating particle memory of %ld MB for type %d\n",npart[type]*sizeof(GENPK_FLOAT_TYPE)/1024/1024,type);
          return -1;
        }
        if (length != get_single_dataset("Masses",masses,length,&hdf_group,fileno)) {
            ret = -1;
            goto exit;
        }
        for(int i = 0; i<npart[type]; i++)
            total_mass_this_file += masses[i];
        *total_mass += total_mass_this_file;
  }
  //Do the final summation here to avoid fp roundoff
  *total_mass += mass[type]*npart[type];
  fieldize(box,field_dims,field,npart[type],pos, masses, mass[type], 1);
exit:
  H5Gclose(hdf_group);
  H5Fclose(hdf_file);
  free(pos);
  return ret;
}
#endif //HDF5

