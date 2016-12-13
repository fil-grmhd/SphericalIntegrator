
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


using namespace SPI;



extern "C" CCTK_INT SphericalIntegrator_RegisterVolumeVar(const CCTK_POINTER_TO_CONST varname,
                                                          const CCTK_POINTER_TO_CONST outname,
                                                          const CCTK_INT sn,
                                                          const CCTK_INT integrate_every)
{
   DECLARE_CCTK_PARAMETERS

   const char* _varname = (char*) varname;
   const char* _outname = (char*) outname;

   assert(_varname != NULL);
   assert(_outname != NULL);

   assert(sn < nslices);
   assert(integrate_every >= 0);

   unsigned num_vol_integrals = 0;
   for(int i = 0; i<slices_1patch.slice().size(); ++i) {
      if(slices_1patch(i,0).integration_type() == volume) {
        num_vol_integrals++;
      }
   }
   if(num_vol_integrals >= max_volume_integrals)
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Too much volume integrals registered, please increase max_volume_integrals parameter to at least %i.",num_vol_integrals+1);

   if(verbose > 1) {
      CCTK_VInfo(CCTK_THORNSTRING,"Registered volume integral on sphere %i:",sn);
      CCTK_VInfo(CCTK_THORNSTRING,"GF: %s | Scalar: %s",_varname,_outname);

      if(integrate_every == 0)
        CCTK_VInfo(CCTK_THORNSTRING,"Never integrating this quantity!");
   }

   // convert varname to lowercase
   // why? CHECK
   string varname_lowercase(_varname);

   for (int i=0; _varname[i]; ++i)
   {
      varname_lowercase[i] = tolower(_varname[i]);
   }

   string outname_lowercase(_outname);

   for (int i=0; _outname[i]; ++i)
   {
      outname_lowercase[i] = tolower(_outname[i]);
   }

   // integration type
   integration_t i_type = volume;
   distrib_method_t d_method = undefined_distrib;

   // return the registered sliced variable number
   return slices_1patch.register_slice(varname_lowercase, outname_lowercase, sn, integrate_every, 0, i_type, d_method, num_vol_integrals);
}


extern "C" CCTK_INT SphericalIntegrator_RegisterSurfaceVar(const CCTK_POINTER_TO_CONST varname,
                                                           const CCTK_POINTER_TO_CONST outname,
                                                           const CCTK_INT sn,
                                                           const CCTK_INT integrate_every,
                                                           const CCTK_INT interpolate_every,
                                                           const CCTK_POINTER_TO_CONST distrib_method)
{
   DECLARE_CCTK_PARAMETERS

   const char* _varname = (char*) varname;
   const char* _outname = (char*) outname;
   const char* _distrib_method = (char*) distrib_method;

   assert(_varname != NULL);
   assert(_distrib_method != NULL);
   assert(_outname != NULL);

   assert(sn < nslices);
   assert(integrate_every >= 0);
   assert(interpolate_every >= 0);

   if(integrate_every < interpolate_every)
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Interpolating less often than integrating for '%s' on sphere %i",_varname,sn);

   if(verbose > 1) {

      CCTK_VInfo(CCTK_THORNSTRING,"Registered surface integral on sphere %i:",sn);
      CCTK_VInfo(CCTK_THORNSTRING,"GF: %s | Scalar: %s",_varname,_outname);

      if(integrate_every == 0)
        CCTK_VInfo(CCTK_THORNSTRING,"Never integrating this quantity!");
      if(interpolate_every == 0)
        CCTK_VInfo(CCTK_THORNSTRING,"Never interpolating this quantity!");
      if((interpolate_every > 0) &&
         (integrate_every > 0) &&
         (interpolate_every != integrate_every))
        CCTK_VInfo(CCTK_THORNSTRING,"CAREFUL: This quantity is interpolated and integrated at different iterations.");
   }


   distrib_method_t d_method = undefined_distrib;
   if (CCTK_Equals(_distrib_method, "const"))
      d_method = constant;
   if (CCTK_Equals(_distrib_method, "single"))
      d_method = single;
   if (CCTK_Equals(_distrib_method, "split"))
      d_method = split;

   if(d_method == undefined_distrib)
      CCTK_WARN(0, "Unkown distribution method! Expecting one defined in spheredata.hh as distrib_t.");

   // convert varname to lowercase
   // why? CHECK
   string varname_lowercase(_varname);

   for (int i=0; _varname[i]; ++i)
   {
      varname_lowercase[i] = tolower(_varname[i]);
   }

   string outname_lowercase(_outname);

   for (int i=0; _outname[i]; ++i)
   {
      outname_lowercase[i] = tolower(_outname[i]);
   }

   // integration type
   integration_t i_type = surface;

   // return the registered sliced variable number
   return slices_1patch.register_slice(varname_lowercase, outname_lowercase, sn, integrate_every, interpolate_every, i_type, d_method, 0);
}
