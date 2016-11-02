
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



extern "C" CCTK_INT SphericalSlice_Register(const CCTK_POINTER_TO_CONST varname, const CCTK_INT sn, const CCTK_INT timelevels, const CCTK_POINTER_TO_CONST distrib_method)
{
   DECLARE_CCTK_PARAMETERS
   
   const char* _varname = (char*) varname;
   const char* _distrib_method = (char*) distrib_method;
   
   assert(_varname != NULL);
   assert(_distrib_method != NULL);
   assert(sn < nslices);
   assert(timelevels > 0);
   
   distrib_method_t d_method;
   if (CCTK_Equals(_distrib_method, "const"))
      d_method = constant;
   if (CCTK_Equals(_distrib_method, "single"))
      d_method = single;
   if (CCTK_Equals(_distrib_method, "split"))
      d_method = split;
   
   if (!CCTK_Equals(_distrib_method, "const") && !CCTK_Equals(_distrib_method, "single") && !CCTK_Equals(_distrib_method, "split"))
      CCTK_WARN(0, "Unkown distribution method!");
   
   // convert varname to lowercase
   string varname_lowercase(_varname);
   
   for (int i=0; _varname[i]; ++i)
   {
      varname_lowercase[i] = tolower(_varname[i]);
   }
   
   if (CCTK_Equals(type[sn], "1patch"))
   {
      // return the registered sliced variable number
      return ONEPATCH_SLICE_IDS+slices_1patch.register_slice(varname_lowercase, sn, timelevels, d_method);
   }
   
   if (CCTK_Equals(type[sn], "2patch"))
   {
      CCTK_WARN(0, "Uh oh....the idea is good but the world isn't ready yet...");
      
      // return the registered sliced variable number
      return TWOPATCH_SLICE_IDS/*+slices_2patch.register_slice(_varname, sn, timelevels)*/;
   }
   
   if (CCTK_Equals(type[sn], "6patch"))
   {
      // return the registered sliced variable number
      return SIXPATCH_SLICE_IDS+slices_6patch.register_slice(varname_lowercase, sn, timelevels, d_method);
   }
   
   
   
   // return -1 on error
   return -1;
}
