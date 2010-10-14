/* Copyright (c) 2010, Simeon Bird <spb41@cam.ac.uk>
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
#define BOOST_TEST_DYN_LINK


#define BOOST_TEST_MODULE MatterPower
#include "gadgetreader.hpp"
#include "gen-pk.h"
#include <math.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#define FLOATS_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( fabs((x) - (y)) <= max<float>(fabs(x),fabs(y))/1e6,(x)<<" is not close to "<<(y))

using namespace GadgetReader;
using namespace std;

BOOST_AUTO_TEST_CASE(check_fieldize)
{
        const int dims=5;
        float field[2*dims*dims*(dims/2+1)]={0};
        float pos[30];
        for(int i=0; i<30; i++)
                pos[i]=i/3.;
        fieldize(10, dims,field,10,10,pos,1);
        //Check off-diagonal elements are zero
        FLOATS_NEAR_TO(field[0],10.7638874);
        BOOST_CHECK_EQUAL(field[3],0);
        BOOST_CHECK_EQUAL(field[20],0);
        BOOST_CHECK_EQUAL(field[125],0);
        //Check on-diagonals
        FLOATS_NEAR_TO(field[124],2.08333);

}

BOOST_AUTO_TEST_CASE(check_invwindow)
{
        FLOATS_NEAR_TO(invwindow(0,3,4,5),71.8177719);
        FLOATS_NEAR_TO(invwindow(4,4,4,5),6111.20801);
        BOOST_CHECK_EQUAL(invwindow(1,1,1,0),0);
}

BOOST_AUTO_TEST_CASE(check_powerspectrum)
{
       float field[2*4*4*(4/2+1)]={1};
       float pow[4];
       int count[4];
       float keffs[4];
       fftwf_complex* outfield;
       outfield=(fftwf_complex *) &field[0];
       fftwf_plan pl=fftwf_plan_dft_r2c_3d(4,4,4,&field[0],outfield, FFTW_ESTIMATE);
       for(int i=0; i<32; i++)
               field[i]+=1;
       BOOST_REQUIRE_EQUAL(powerspectrum(4,&pl,outfield,4,pow,count,keffs),0);
       FLOATS_NEAR_TO(keffs[2],2.43421054);
       BOOST_CHECK_EQUAL(count[1],26);
       BOOST_CHECK_EQUAL(count[0],1);
       FLOATS_NEAR_TO(pow[0],0.129150391);
       FLOATS_NEAR_TO(pow[1],0.0134053277);
       FLOATS_NEAR_TO(pow[2],0.0139268516);
       FLOATS_NEAR_TO(pow[3],0.035415493);
       fftwf_destroy_plan(pl);
}

BOOST_AUTO_TEST_CASE(check_read_fieldize)
{
        float field[2*4*4*(4/2+1)];
        GSnap snap("test_g2_snap", false);
        BOOST_REQUIRE_MESSAGE(snap.GetNumFiles()==2,"Failed to find test snapshot data");
        BOOST_REQUIRE_EQUAL(read_fieldize(field, &snap, BARYON_TYPE, 3000,4),0);
        FLOATS_NEAR_TO(field[10],0);
        FLOATS_NEAR_TO(field[0],1.64618);
        FLOATS_NEAR_TO(field[15],0.700735);
        BOOST_REQUIRE_EQUAL(read_fieldize(field, &snap,BULGE_TYPE, 3000,4),1);
}

BOOST_AUTO_TEST_CASE(check_type_str)
{
        BOOST_CHECK_EQUAL(type_str(BARYON_TYPE),"by");
        BOOST_CHECK_EQUAL(type_str(DM_TYPE),"DM");
}

BOOST_AUTO_TEST_CASE(check_nexttwo)
{
        BOOST_CHECK_EQUAL(nexttwo(12),16);
        BOOST_CHECK_EQUAL(nexttwo(8),8);
        BOOST_CHECK_EQUAL(nexttwo(1),1);
}
