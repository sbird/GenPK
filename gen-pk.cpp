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
 * Support "fake" kspace neutrinos.
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
#include <math.h>
//For getopt
#include <unistd.h>
//For omp_get_num_procs
#include <omp.h>
#include <stdlib.h>
#include <sstream>

/** Maximal size of FFT grid. 
 * In practice 1024 means we need just over 4GB, as sizeof(float)=4*/
#define FIELD_DIMS 1024

using namespace GadgetReader;
using namespace std;

#ifndef NOHDF5
#include <hdf5.h>

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
     if (H5Fis_hdf5(fname.c_str()) <= 0)
        fname = std::string();
  }

  /* Restore previous error handler */
  H5Eset_auto(error_stack, old_func, old_client_data);
  return fname;
}
#endif

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

/** \file 
 * File containing main() */

/** Main function. Accepts arguments, uses GadgetReader to open the snapshot, prints some header info, 
 * allocates FFT memory and creates a plan, calls read_fieldize and powerspectrum, prints the P(k) 
 * to a file, and then frees the memory*/
int main(int argc, char* argv[]){
  int nrbins,field_dims=0,type;
  float *field, *power, *keffs;
  int *count;
  int64_t npart_total[N_TYPE];
  string infiles(""),outdir("");
  char c;
  double box;
  GSnap * snap = NULL;
  int fileno=0;
  bool use_hdf5 = false;
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
  unsigned i_fileno=0;
#ifndef NOHDF5
    /*ffname is a copy of input filename for extension*/
    /*First open first file to get header properties*/
    std::string fname = find_first_hdf_file(infiles);
    std::string ffname = fname;
    if ( !fname.empty() ){
            /*See if we have been handed the first file of a set:
             * our method for dealing with this closely mirrors
             * HDF5s family mode, but we cannot use this, because
             * our files may not all be the same size.*/
	    i_fileno = fname.find(".0.hdf5")+1;
        use_hdf5 = true;
        double atime, redshift, h100;
        //Get the header and print out some useful things
        if(load_hdf5_header(fname.c_str(), &atime, &redshift, &box, &h100, npart_total)) {
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
      for(type=0;type<N_TYPE;type++)
          npart_total[type]=snap->GetNpart(type);
      //Get the header and print out some useful things
    box=snap->GetHeader().BoxSize;
    fprintf(stderr, "Boxsize=%g, ",box);
    fprintf(stderr, "NPart=(%g,%g,%g,%g,%g,%g)**3\n",cbrt(npart_total[0]),cbrt(npart_total[1]),cbrt(npart_total[2]),cbrt(npart_total[3]),cbrt(npart_total[4]),cbrt(npart_total[5]));
    fprintf(stderr, "Masses=[%g %g %g ]\n",snap->GetHeader().mass[0],snap->GetHeader().mass[1],snap->GetHeader().mass[2]);
    fprintf(stderr, "redshift=%g, Ω_M=%g\n",snap->GetHeader().redshift,snap->GetHeader().Omega0);
  }
  //Work out how large a field we need
  for(type=0;type<N_TYPE;type++){
    int tmp=2*nexttwo(cbrt(npart_total[type]));
    field_dims=std::max(field_dims, std::min(tmp, FIELD_DIMS));
  }
  nrbins=floor(sqrt(3)*((field_dims+1.0)/2.0)+1);
  //Memory for the field
  /* Allocating a bit more memory allows us to do in-place transforms.*/
  if(!(field=(float *)fftwf_malloc(2*field_dims*field_dims*(field_dims/2+1)*sizeof(float)))){
  	fprintf(stderr,"Error allocating memory for field\n");
  	return 1;
  }
  string filename=outdir;
  size_t last=infiles.find_last_of("/\\");
  //Set up FFTW
  outfield=(fftwf_complex *) &field[0];
  if(!fftwf_init_threads()){
  		  fprintf(stderr,"Error initialising fftw threads\n");
  		  return 0;
  }
  fftwf_plan_with_nthreads(omp_get_num_procs());
  pl=fftwf_plan_dft_r2c_3d(field_dims,field_dims,field_dims,&field[0],outfield, FFTW_ESTIMATE);
  //Allocate memory for output
  power=(float *) malloc(nrbins*sizeof(float));
  count=(int *) malloc(nrbins*sizeof(int));
  keffs=(float *) malloc(nrbins*sizeof(float));
  if(!power || !count || !keffs){
  	fprintf(stderr,"Error allocating memory for power spectrum.\n");
        return 1;
  }
  /*Now make a power spectrum for each particle type*/
  for(type=0; type<N_TYPE; type++){
        if(npart_total[type] == 0)
            continue;
        if (use_hdf5){
            //Find filename
            if(i_fileno != std::string::npos){
		        std::ostringstream convert;
		        convert<<fileno;
                ffname = fname.replace(i_fileno, 1, convert.str());
		    }
            else
             break;
            /*If we ran out of files, we're done*/
            if(!(file_readable(ffname.c_str()) && H5Fis_hdf5(ffname.c_str()) > 0))
                    break;
            read_fieldize_hdf5(field, ffname.c_str(), type, box, field_dims, fileno);
            fileno++;
        }
        else{
            if(read_fieldize(field,snap,type, box, field_dims))
                continue;
        }
        if(powerspectrum(field_dims,&pl,outfield,nrbins, power,count,keffs))
                continue;
        filename=outdir;
        filename+="/PK-"+type_str(type)+"-"+infiles.substr(last+1);
        print_pk(filename,nrbins,keffs,power,count);
  }

  //Free memory
  delete snap;
  free(power);
  free(count);
  free(keffs);
  fftwf_free(field);
  fftwf_destroy_plan(pl);
  return 0;
}

