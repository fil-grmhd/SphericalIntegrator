
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



extern "C" CCTK_REAL SphericalIntegrator_SurfaceIntegrate(const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT varno)
{
   DECLARE_CCTK_PARAMETERS

   assert(varno >= 0);

   assert(varno < slices_1patch.slice().size());

   cGH const * restrict const cctkGH = (cGH*) cctkGH_;

   if(verbose > 0) {
     CCTK_VInfo(CCTK_THORNSTRING,"Taking manual surface integral of '%s' on sphere %i.",slices_1patch(varno,0).varname().c_str(),slices_1patch(varno,0).ID());
   }

   // return the surface integral
   return slices_1patch(varno, 0).integrate_surface(cctkGH);
}

extern "C" CCTK_REAL SphericalIntegrator_VolumeIntegrate(const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT varno)
{
   DECLARE_CCTK_PARAMETERS

   assert(varno >= 0);

   cGH const * restrict const cctkGH = (cGH*) cctkGH_;

   if(verbose > 0) {
     CCTK_VInfo(CCTK_THORNSTRING,"Taking manual volume integral of '%s' over sphere %i.",slices_1patch(varno,0).varname().c_str(),slices_1patch(varno,0).ID());
   }

   const CCTK_INT sum_reduction_handle = CCTK_ReductionHandle("sum");
   if(sum_reduction_handle < 0)
       CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"SphericalIntegrator: Unable to get reduction handel 'sum'.");

   assert(varno < slices_1patch.slice().size());

   // return the surface integral
   return slices_1patch(varno, 0).integrate_volume(cctkGH,sum_reduction_handle);
}

extern "C" void SphericalIntegrator_CollectiveIntegration(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose > 0) {
    CCTK_INFO("Computing integrals.");
  }

  // get reduction handle for volume integration
  const CCTK_INT sum_reduction_handle = CCTK_ReductionHandle("sum");
  if(sum_reduction_handle < 0)
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"SphericalIntegrator: Unable to get reduction handel 'sum'.");
  // precompute delta volume
  CCTK_REAL dx3 = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];

  vector<CCTK_INT> vol_vars;
  // collect volume integration vars, to get tmp gf indices
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    if(slices_1patch(i,0).integration_type() == volume) {
      if(slices_1patch(i,0).has_constant_radius())
        vol_vars.push_back(i);
      else
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Volume integrals over elliptical surfaces not supported yet. (sn = %i)",slices_1patch(i,0).ID());
    }
  }

  // get GF index size to shift indices in GF arrays
  const CCTK_INT gf_size = UTILS_GFSIZE(cctkGH);

  // loop over all vars
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    // integrate only if it is necessary
    if(cctk_iteration % slices_1patch(i,0).integrate_every() == 0) {
      // integrate and store result
      if(verbose > 0) {
        CCTK_VInfo(CCTK_THORNSTRING,"Integrating '%s' in collective mode on sphere %i.",slices_1patch(i,0).varname().c_str(),slices_1patch(i,0).ID());
      }

      if(slices_1patch(i,0).integration_type() == surface)
        slices_1patch(i, 0).integrate_surface(cctkGH);
      else {
        slices_1patch(i,0).integrate_volume(cctkGH,sum_reduction_handle);
      }
    }
  }



}
