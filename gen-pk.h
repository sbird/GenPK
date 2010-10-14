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

#ifndef GEN_PK_HEADER
#define GEN_PK_HEADER

#include <fftw3.h>
#include "gadgetreader.hpp"
#include <string>

int nexttwo(int);
int print_pk(std::string filename, int nrbins, float * keffs, float * power, int* count);
void help(void);
std::string type_str(int type);

#ifdef __cplusplus
extern "C" {
#endif
/* Fieldize. positions should be an array of size 3*particles 
 * (like the output of read_gadget_float3)
 * out is an array of size [dims*dims*dims]*/
/* the "extra" switch, if set to one, will assume that the output 
 * is about to be handed to an FFTW in-place routine, 
 * and make it skip the last 2 places of the each row in the last dimension */
int fieldize(double boxsize, int dims, float *out, int total_particles, int segment_particles, float *positions, int extra);

/* The window function associated with the above.*/
float invwindow(int kx, int ky, int kz, int n);
int read_fieldize(float * field, GadgetReader::GSnap* snap, int type, double box, int field_dims);
/*Note we will need some contiguous memory space after the actual data in field.
 * The real input data has size
 * dims*dims*dims
 * The output has size dims*dims*(dims/2+1) *complex* values
 * So need 2*dims*dims*(dims/2+1) float space.
 * Need at least floor(sqrt(3)*abs((dims+1.0)/2.0)+1) values in power and count.*/
int powerspectrum(int dims, fftwf_plan* pl,fftwf_complex* outfield, int nrbins, float *power, int *count, float *keffs);
#ifdef __cplusplus
}
#endif

#endif
