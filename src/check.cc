
/* Copyright 2013 Peter Diener, Nils Dorband, Roland Haas, Ian Hinder,
Christian Ott, Denis Pollney, Thomas Radke, Christian Reisswig, Erik
Schnetter, Barry Wardell and Burkhard Zink

This file is part of Llama.

Llama is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License, or (at your
option) any later version.

Llama is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Llama.  If not, see <http://www.gnu.org/licenses/>. */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "slices.hh"



namespace SPI {



using namespace SPI;


extern "C" void SphericalIntegrator_CheckAndUpdate(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   if(verbose>1)
     CCTK_VInfo(CCTK_THORNSTRING,"Checking and updating slices (it=%i).",cctk_iteration);

   for (int n=0; n < nslices; ++n)
   {
      if (ss_valid[n] > 0 && ss_active[n] == 0)
      {
        CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Slice #%d has ss_valid set to a positive value, but does "
                     "not have ss_active set.  This is an error in the thorn "
                     "which calculated this surface", n);
      }
   }
   // update origins of spheres if bnstracker is used
   if(bnstracker_positions) {
      // try to get bnstracker var indices, check if they are there
      CCTK_INT star1_index_x = CCTK_VarIndex("BNSTrackerGen::bns_x_1");
      CCTK_INT star1_index_y = CCTK_VarIndex("BNSTrackerGen::bns_y_1");

      CCTK_INT star2_index_x = CCTK_VarIndex("BNSTrackerGen::bns_x_2");
      CCTK_INT star2_index_y = CCTK_VarIndex("BNSTrackerGen::bns_y_2");

      CCTK_INT star1_index_x_md = CCTK_VarIndex("BNSTrackerGen::bns_x_md_1");
      CCTK_INT star1_index_y_md = CCTK_VarIndex("BNSTrackerGen::bns_y_md_1");

      CCTK_INT star2_index_x_md = CCTK_VarIndex("BNSTrackerGen::bns_x_md_2");
      CCTK_INT star2_index_y_md = CCTK_VarIndex("BNSTrackerGen::bns_y_md_2");

      // checking one index should be enough
      if(star1_index_x < 0)
        CCTK_WARN(0,"Can't access BNSTrackerGen variables: Is it activated?");

      // get the pointer to the bnstracker vars
      CCTK_REAL* star1_ptr_x = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star1_index_x,NULL);
      CCTK_REAL* star1_ptr_y = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star1_index_y,NULL);

      CCTK_REAL* star2_ptr_x = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star2_index_x,NULL);
      CCTK_REAL* star2_ptr_y = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star2_index_y,NULL);

      CCTK_REAL* star1_ptr_x_md = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star1_index_x_md,NULL);
      CCTK_REAL* star1_ptr_y_md = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star1_index_y_md,NULL);

      CCTK_REAL* star2_ptr_x_md = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star2_index_x_md,NULL);
      CCTK_REAL* star2_ptr_y_md = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,star2_index_y_md,NULL);

      for(int i = 0; i<nslices; ++i) {
         switch(ss_track[i]) {
            case 0: {
              // do nothing, no update wanted
              break;
            }
            case 1: {
              // track star 1
              ss_origin_x[i] = *star1_ptr_x;
              ss_origin_y[i] = *star1_ptr_y;
              ss_origin_z[i] = 0;
              break;
            }
            case 2: {
              // track star 2
              ss_origin_x[i] = *star2_ptr_x;
              ss_origin_y[i] = *star2_ptr_y;
              ss_origin_z[i] = 0;
              break;
            }
            case 3: {
              // track star 1 md
              ss_origin_x[i] = *star1_ptr_x_md;
              ss_origin_y[i] = *star1_ptr_y_md;
              ss_origin_z[i] = 0;
              break;
            }
            case 4: {
              // track star 2 md
              ss_origin_x[i] = *star2_ptr_x_md;
              ss_origin_y[i] = *star2_ptr_y_md;
              ss_origin_z[i] = 0;
              break;
            }
         }
      }
   }
   // go through all 1-patch slices and check and change most recent timelevel
   for (int j=0; j < slices_1patch.slice().size(); ++j)
   {
      int id = slices_1patch.slice()[j].front().ID(); // n-th slice in par-file

      vect<CCTK_REAL,3> new_origin(ss_origin_x[id], ss_origin_y[id], ss_origin_z[id]);
      slices_1patch.slice()[j].front().origin() = new_origin;
   }

   // TODO: also set other ss_info variables!
   //       But on the other hand calculating the surface area at each timestep might be expensive....
}

}
