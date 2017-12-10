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

/** \mainpage
 * \section intro_sec Introduction
 * MatPow is a 3D matter power spectrum estimator, written in C++ and parallelised with OpenMP.
 *
 * \section feat_sec Features
 *
 * A matter power spectrum estimator, which can read most Gadget-I and II format files through the use of the
 * GadgetReader library. Calculates the matter power spectrum quickly, as well as the expected error from
 * cosmic variance.
 *
 * Outputs one file per particle type found.
 *
 * \section missing_sec Missing Features
 *
 * Should rebin to ensure enough modes per bin
 *
 * Use Volker Springel's fold-on-itself to calculate small-scale power below the Nyquist frequency
 *
 * \section usage_sec Usage
 *
 * See:
 * ./gen-pk --help
 *
 * \section req_sec Requirements
 * A C++ compiler with map, vector, set getopt, and stdint.h
 *
 * GadgetReader library, preferably installed in ../GadgetReader
 * (edit the first line of the Makefile for other locations)
 *
 * FFTW3
 *
 * Boost::Test library >= 1.34 for the test suite
 */
#include "gen-pk.h"
#include <string.h>
#include <math.h>
//For getopt
#include <unistd.h>
//For omp_get_num_procs
#ifdef _OPENMP
#include <omp.h>
#else
#warning "You are compiling with OpenMP disabled (probably Apple clang). This will be slow! Try installing gcc."
#endif
#include <stdlib.h>

/** Maximal size of FFT grid.
 * In practice 1024 means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 3072L

using namespace GadgetReader;
using namespace std;


/** \file
 * File containing main() */

/** Main function. Accepts arguments, uses GadgetReader to open the snapshot, prints some header info,
 * allocates FFT memory and creates a plan, calls read_fieldize and powerspectrum, prints the P(k)
 * to a file, and then frees the memory*/
int main(int argc, char* argv[])
{
  int64_t field_dims=0;
  int type;
  GENFLOAT *field;
  double *power, *keffs;
  int *count;
  int64_t npart_total[N_TYPE];
  double mass[N_TYPE];
  string infiles(""),jinfiles(""),outdir("");
  char c;
  int crosstype = -1;
  bool stars_are_baryons = false;
  double box;
  double Omega0;
  GSnap * snap = NULL;
  bool use_hdf5 = false;
  bool use_bigfile = false;
  fftw_plan pl;
  fftw_complex *outfield;
  while((c = getopt(argc, argv, "i:j:o:c:s:h")) !=-1){
    switch(c){
        case 'o':
           outdir=static_cast<string>(optarg);
           break;
        case 'i':
           infiles=static_cast<string>(optarg);
           break;
        case 'j':
	   jinfiles=static_cast<string>(optarg);
           break;
        case 'c':
           crosstype = static_cast<int>(atoi(optarg));
           break;
        case 's':
            stars_are_baryons = static_cast<bool>(atoi(optarg));
            break;
        case 'h':
        default:
           help();
           return 0;
      }
  }
#ifndef NOHDF5
  /*ffname is a copy of input filename for extension*/
  /*First open first file to get header properties*/
  std::vector<std::string> fnames = find_hdf_set(infiles);
  if ( !fnames.empty() ){
      fprintf(stderr, "Found %lu hdf5 files in snapshot\n",fnames.size());
      /*See if we have been handed the first file of a set:
       * our method for dealing with this closely mirrors
       * HDF5s family mode, but we cannot use this, because
       * our files may not all be the same size.*/
      use_hdf5 = true;
      double atime, redshift, h100;
      //Get the header and print out some useful things
      if(load_hdf5_header(fnames[0].c_str(), &atime, &redshift, &box, &h100, npart_total, mass)) {
        fprintf(stderr, "Could not load header\n");
        return 1;
      }
  }
  else
#endif
#ifndef NOBIGFILE
  if(is_bigfile(infiles.c_str())) {
      use_bigfile = true;
      double atime, redshift, h100;
      //Get the header and print out some useful things
      if(load_bigfile_header(infiles.c_str(), &atime, &redshift, &box, &h100, npart_total, mass, &Omega0)) {
        fprintf(stderr, "Could not load header\n");
        return 1;
      }
  }
  else
#endif
  {
    //Open the snapshot
    snap = new GSnap(infiles);
    if(outdir.empty() || (!use_hdf5 && snap->GetNumFiles() < 1)){
            help();
            return 0;
    }
      for(type=0;type<N_TYPE;type++){
          npart_total[type]=snap->GetNpart(type);
          mass[type] = snap->GetHeader().mass[type];
      }
      //Get the header and print out some useful things
    box=snap->GetHeader().BoxSize;
    printf("Boxsize=%g, ",box);
    printf("NPart=(%g,%g,%g,%g,%g,%g)**3\n",cbrt(npart_total[0]),cbrt(npart_total[1]),cbrt(npart_total[2]),cbrt(npart_total[3]),cbrt(npart_total[4]),cbrt(npart_total[5]));
    printf("Masses=[%g %g %g ]\n",snap->GetHeader().mass[0],snap->GetHeader().mass[1],snap->GetHeader().mass[2]);
    printf("redshift=%g, Ω_M=%g\n",snap->GetHeader().redshift,snap->GetHeader().Omega0);
  }
  //Work out how large a field we need
  for(type=0;type<N_TYPE;type++){
    int64_t tmp=2*nexttwo(cbrt(npart_total[type]));
    field_dims=std::max(field_dims, std::min(tmp, FIELD_DIMS));
  }
  const int nrbins=field_dims;
  //Memory for the field
  printf("FFT grid dimension: %lu\n",field_dims);
  const size_t field_size = 2*field_dims*(size_t) field_dims*(field_dims/2+1);
  /* Allocating a bit more memory allows us to do in-place transforms.*/
  if(!(field=(GENFLOAT *)fftw_malloc(field_size*sizeof(GENFLOAT)))){
  	fprintf(stderr,"Error allocating memory for field\n");
  	return 1;
  }
  string filename=outdir;
  size_t last=infiles.find_last_of("/\\");
  //Set up FFTW
  outfield=(fftw_complex *) &field[0];
  if(!fftw_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  return 0;
  }
#ifdef _OPENMP
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  pl=fftw_plan_dft_r2c_3d(field_dims,field_dims,field_dims,&field[0],outfield, FFTW_ESTIMATE);
  //Allocate memory for output
  power=(double *) malloc(nrbins*sizeof(double));
  count=(int *) malloc(nrbins*sizeof(int));
  keffs=(double *) malloc(nrbins*sizeof(double));
  if(!power || !count || !keffs){
  	fprintf(stderr,"Error allocating memory for power spectrum.\n");
        return 1;
  }
  if  ((crosstype < 0) && (jinfiles.size()==0)) {
    /*Now make a power spectrum for each particle type*/
    for(type=0; type<N_TYPE; type++){
          //Zero field
          if(npart_total[type] == 0)
              continue;
          memset(field, 0, field_size*sizeof(GENFLOAT));
          double total_mass = 0;
          if (use_hdf5){
              for(unsigned fileno = 0; fileno < fnames.size(); ++fileno)
                  /*If we couldn't read the coordinate array for some reason, go on to the next type*/
                  if(read_fieldize_hdf5(field, fnames[fileno].c_str(), type, box, field_dims, &total_mass, fileno) < 0)
                      continue;
              if(type == 0 && stars_are_baryons) {
                for(unsigned fileno = 0; fileno < fnames.size(); ++fileno)
                  read_fieldize_hdf5(field, fnames[fileno].c_str(), STARS_TYPE, box, field_dims, &total_mass, fileno);
              }
          }
          else if(use_bigfile) {
                  read_fieldize_bigfile(field, infiles.c_str(), type, box, field_dims, &total_mass, npart_total, mass, Omega0);
                  if(type == 0 && stars_are_baryons) {
                      read_fieldize_bigfile(field, infiles.c_str(), STARS_TYPE, box, field_dims, &total_mass, npart_total, mass, Omega0);
                  }
          }
          else {
              read_fieldize(field,snap,type, box, field_dims, &total_mass);
              if(type == 0 && stars_are_baryons) {
                  read_fieldize(field,snap,STARS_TYPE, box, field_dims, &total_mass);
              }
          }
          printf("total_mass in type %d = %g\n", type, total_mass);
          fftw_execute(pl);
          if(powerspectrum(field_dims,outfield, outfield, nrbins, power,count,keffs, total_mass, total_mass))
                  continue;
          filename=outdir;
          filename+="/PK-"+type_str(type)+"-"+infiles.substr(last+1);
          print_pk(filename,nrbins,keffs,power,count);
    }
  } else if (jinfiles.size()>0)  {
    // do cross correlation across two files
    bool use_hdf52=false;
    GSnap * snap2 = NULL;
    #ifndef NOHDF5
    std::vector<std::string> fnames2 = find_hdf_set(jinfiles);
    if ( !fnames2.empty() ){
        /*See if we have been handed the first file of a set:
         * our method for dealing with this closely mirrors
         * HDF5s family mode, but we cannot use this, because
         * our files may not all be the same size.*/
        use_hdf52 = true;
    }
    else
    #endif
        snap2 = new GSnap(jinfiles);
    /*Now make a power spectrum for each particle type*/
    for(type=0; type<N_TYPE; type++){
          if(npart_total[type] == 0)
              continue;
          //Memory for the field
          /* Allocating a bit more memory allows us to do in-place transforms.*/
          GENFLOAT * field2;
          if(!(field2=(GENFLOAT *)fftw_malloc(field_size*sizeof(GENFLOAT)))){
            fprintf(stderr,"Error allocating memory for second field\n");
            return 1;
          }
          memset(field, 0, field_size*sizeof(GENFLOAT));
          memset(field2, 0, field_size*sizeof(GENFLOAT));
          fftw_complex * outfield2=(fftw_complex *) &field2[0];
          fftw_plan pl2=fftw_plan_dft_r2c_3d(field_dims,field_dims,field_dims,&field2[0],outfield2, FFTW_ESTIMATE);

          fprintf(stderr,"Reading...\n");
          //Get the new filename if HDF5
          double total_mass = 0;
          if (use_hdf52){
              for(unsigned fileno = 0; fileno < fnames.size(); ++fileno)
                  read_fieldize_hdf5(field, fnames[fileno].c_str(), 1, box, field_dims, &total_mass, fileno);
          }
          else if(use_bigfile) {
                  read_fieldize_bigfile(field, infiles.c_str(), type, box, field_dims, &total_mass, npart_total, mass, Omega0);
          }
          else if(read_fieldize(field,snap,type, box, field_dims, &total_mass))
            continue;
          double total_mass2 = 0;
          //Get the new filename if HDF5
          if (use_hdf52){
              for(unsigned fileno = 0; fileno < fnames2.size(); ++fileno)
                  read_fieldize_hdf5(field2, fnames2[fileno].c_str(), 1, box, field_dims, &total_mass2, fileno);
          }
          else if(use_bigfile) {
                  read_fieldize_bigfile(field, infiles.c_str(), type, box, field_dims, &total_mass, npart_total, mass, Omega0);
          }
          else if(read_fieldize(field2,snap2,type, box, field_dims, &total_mass2))
            continue;
          fftw_execute(pl);
          fftw_execute(pl2);
          if(powerspectrum(field_dims,outfield, outfield2, nrbins, power,count,keffs, total_mass, total_mass2))
                  continue;
          filename=outdir;
          filename+="/PX-"+type_str(type)+"-"+infiles.substr(last+1);
          print_pk(filename,nrbins,keffs,power,count);
    }
  }
  else
  {
      //Do a cross-correlation inside the file
        if(crosstype > N_TYPE || npart_total[1] == 0 || npart_total[crosstype] == 0){
            fprintf(stderr, "Can't cross-correlate types not present in snapshot\n");
            return 1;
        }
       //Memory for the field
       /* Allocating a bit more memory allows us to do in-place transforms.*/
       GENFLOAT * field2;
       if(!(field2=(GENFLOAT *)fftw_malloc(field_size*sizeof(GENFLOAT)))){
       	fprintf(stderr,"Error allocating memory for second field\n");
       	return 1;
       }
       memset(field, 0, field_size*sizeof(GENFLOAT));
       memset(field2, 0, field_size*sizeof(GENFLOAT));
       fftw_complex * outfield2=(fftw_complex *) &field2[0];
       fftw_plan pl2=fftw_plan_dft_r2c_3d(field_dims,field_dims,field_dims,&field2[0],outfield2, FFTW_ESTIMATE);
       //Get the DM
       double total_mass = 0;
       if (use_hdf5){
           for(unsigned fileno = 0; fileno < fnames.size(); ++fileno)
               read_fieldize_hdf5(field, fnames[fileno].c_str(), 1, box, field_dims, &total_mass, fileno);
       }
       else if(use_bigfile) {
               read_fieldize_bigfile(field, infiles.c_str(), 1, box, field_dims, &total_mass, npart_total, mass, Omega0);
       }
       else
           read_fieldize(field,snap,1, box, field_dims, &total_mass);
       //Get the other species
       double total_mass2 = 0;
       if (use_hdf5){
           for(unsigned fileno = 0; fileno < fnames.size(); ++fileno)
               read_fieldize_hdf5(field2, fnames[fileno].c_str(), crosstype, box, field_dims, &total_mass2, fileno);
       }
       else if(use_bigfile) {
               read_fieldize_bigfile(field2, infiles.c_str(), crosstype, box, field_dims, &total_mass2, npart_total, mass, Omega0);
       }
       else
           read_fieldize(field2,snap,crosstype, box, field_dims, &total_mass2);
       //Do FFT of DM
	   fftw_execute(pl);
       //Do FFT of other species
	   fftw_execute(pl2);
       powerspectrum(field_dims,outfield, outfield2, nrbins, power,count,keffs, total_mass, total_mass2);
       filename=outdir;
       filename+="/PK-DMx"+type_str(crosstype)+"-"+infiles.substr(last+1);
       //Multiply by a mass factor
//        for (int i=0; i< nrbins; i++) {
//            power[i]*= mass[1]/mass[crosstype];
//        }
       print_pk(filename,nrbins,keffs,power,count);
  }

  //Free memory
  delete snap;
  free(power);
  free(count);
  free(keffs);
  fftw_free(field);
  fftw_destroy_plan(pl);
  return 0;
}

