/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
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

/** \file
 * GenPK Header*/
#ifndef GEN_PK_HEADER
#define GEN_PK_HEADER

#include <fftw3.h>
#include "gadgetreader.hpp"
#include <string>

#ifdef DOUBLE_PRECISION
    #define FLOAT_TYPE double
#else
    #define FLOAT_TYPE float
#endif

/** Returns the next power of two greater than its argument*/
int nexttwo(int);

/** Function to print the calculated power spectrum to a file
 * Will print three columns, k_eff P(k) and N_modes.
 * @param filename File to print to
 * @param nrbins  Number of bins to print
 * @param keffs Effective k; this is the center of the bins, weighted by the fact that there are more 
 * modes at higher k
 * @param power Power spectrum
 * @param count Number of modes in each bin*/
int print_pk(std::string filename, int nrbins, double * keffs, double * power, int* count);

/** Prints usage info */
void help(void);

/** Function to convert particle type to a string
 * @param type Particle type, such as BARYON_TYPE*/
std::string type_str(int type);


/** Wrapper function to read the data and pass it to fieldize in chunks.
 * Returns 1 is some error occured allocating memory, and zero if successful.
 * @param field Pointer to memory for output field
 * @param snap GadgetReader object for snapshot
 * @param type Particle type to read; since POS always has all types, complications are minimal
 * @param box Boxsize, in kpc
 * @param field_dims Length of side of the field above. */
int read_fieldize(float * field, GadgetReader::GSnap* snap, int type, double box, int field_dims, double * total_mass);

/** Wrapper function to read the data from an HDF5 snapshot and pass it to fieldize.
 * Returns 1 is some error occured allocating memory, and zero if successful.
 * @param field Pointer to memory for output field
 * @param ffname filename of hdf snapshot
 * @param type Particle type to read
 * @param box Boxsize, in kpc
 * @param field_dims Length of side of the field above.
 * @param fileno number of file to use.
 * @param total_mass Variable to store total mass of particles */
int read_fieldize_hdf5(float * field, const char *ffname, int type, double box, int field_dims, double * total_mass,int fileno);

/* this routine loads header data from the first file of an HDF5 snapshot.*/
int load_hdf5_header(const char *ffname, double  *atime, double *redshift, double *box100, double *h100, int64_t *npart_all, double * mass);

/** Finds a snapshot set of hdf files from the given initial filename
 */
std::vector<std::string> find_hdf_set(const std::string& infname);

/** Takes an array of particle positions and places them into a grid, ready for use by FFTW.
 * OpenMP-parallel. Returns 0. Uses a Cloud-In-Cell algorithm.
 * @param boxsize The size of the box, in kpc.
 * @param dims The size of the grid to use
 * @param out Pointer to output array, of size dims**3. Should be already initialised. 
 * @param total_particles Total number of particles that will be placed on the grid. Used for density calculation.
 * @param segment_particles Number of particles to place on the grid in this call. 
 * @param positions Array of particle positions, of size 3*segment_particles (like the output of GetBlocks)
 * @param masses Array containing particle masses, if variable. Null if particle masses are constant.
 * @param mass Particle mass if constant for all of this type.
 * @param extra If this is 1, assume that the output is about to be handed to an FFTW in-place routine, 
 * and make it skip the last 2 places of the each row in the last dimension */
int fieldize(double boxsize, int dims, float *out, int64_t segment_particles, FLOAT_TYPE *positions, FLOAT_TYPE * masses, double mass, int extra);

#ifdef __cplusplus
extern "C" {
#endif
/** The inverse window function associated with the Cloud-in-cell algorithm.
 * The power spectrum needs to deconvolve this.
 * @param kx k_x
 * @param ky k_y
 * @param kz k_z
 * @param n Total size of grid */
float invwindow(int kx, int ky, int kz, int n);
/** Wrapper around FFTW to calculate the 3D power spectrum from a 3D field.
 * Returns 0.
 * @param dims Size of grid to FFT.
 * @param outfield Pointer to memory where the FFT stores its output, as was specified to the FFTW plan.
 * This is likely to be an in-place transform, in which case 
 * we will need some contiguous memory space after the actual data in field.
 * The real input data has size dims**3
 * The output has size dims*dims*(dims/2+1) *complex* values
 * So we must allocate 2*dims*dims*(dims/2+1)*sizeof(float).
 * @param outfield2 Pointer to second array, just like outfield. What is computed is outfield*outfield2. May be the same array.
 * @param nrbins Number of bins in the output.
 * @param power Pointer to memory to output powerspectrum to. Needs to store nrbins values.
 * @param count Ditto for modes per bin
 * @param keffs Ditto for effective k.*/
int powerspectrum(int dims, fftwf_complex* outfield, fftwf_complex* outfield2, int nrbins, double *power, int *count, double *keffs, double total_mass, double total_mass2);

#ifdef __cplusplus
}
#endif

#endif
