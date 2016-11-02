
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


using namespace SPS;



extern "C" CCTK_REAL SphericalSlice_Integrate(const CCTK_INT varno, const CCTK_INT timelevel)
{
   DECLARE_CCTK_PARAMETERS
   
   assert(timelevel >= 0);
   assert(varno >= 0);
   
   CCTK_REAL result = 0;
   
   if (is_1patch(varno))
   {
      assert(INDEX1P(varno) < slices_1patch.slice().size());
      result = slices_1patch(INDEX1P(varno), timelevel).integrate();
   }
   
   if (is_2patch(varno))
      CCTK_WARN(0, "Uh oh....the idea is good but the world isn't ready yet...");
   
   if (is_6patch(varno))
   {
      assert(INDEX6P(varno) < slices_6patch.slice().size());
      result = slices_6patch(INDEX6P(varno), timelevel).integrate();
   }
   
   // return the surface integral
   return result;
}
