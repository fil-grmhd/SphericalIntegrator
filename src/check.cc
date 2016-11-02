
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



namespace SPS {



using namespace SPS;


extern "C" void SphericalSlice_CheckAndUpdate(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
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
   
   // go through all 1-patch slices and check and change most recent timelevel
   for (int j=0; j < slices_1patch.slice().size(); ++j)
   {
      int id = slices_1patch.slice()[j].front().ID(); // n-th slice in par-file
      
      vect<CCTK_REAL,3> new_origin(ss_origin_x[id], ss_origin_y[id], ss_origin_z[id]);
      slices_1patch.slice()[j].front().origin() = new_origin;
   }
   
   
   // go through all 6-patch slices and check and change most recent timelevel
   for (int j=0; j < slices_6patch.slice().size(); ++j)
   {
      int id = slices_6patch.slice()[j].front().ID(); // n-th slice in par-file
      
      vect<CCTK_REAL,3> new_origin(ss_origin_x[id], ss_origin_y[id], ss_origin_z[id]);
      slices_6patch.slice()[j].front().origin() = new_origin;
   }
   
   // TODO: also set other ss_info variables!
   //       But on the other hand calculating the surface area at each timestep might be expensive....
   
}


}