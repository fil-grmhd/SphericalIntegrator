
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



extern "C" CCTK_REAL SphericalIntegrator_Integrate(const CCTK_INT varno, const CCTK_INT timelevel)
{
   DECLARE_CCTK_PARAMETERS

   assert(timelevel >= 0);
   assert(varno >= 0);

   CCTK_REAL result = 0;

   assert(varno < slices_1patch.slice().size());
   // return the surface integral
   return slices_1patch(varno, timelevel).integrate();
}

extern "C" void SphericalIntegrator_CollectiveIntegration(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  vector<CCTK_INT> vars[nslices];

  // collect all vars and sort them to their sphere ID
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    // only sync if we also want to integrate
    if(cctk_iteration % slices_1patch(i,0).integrate_every() == 0)
      vars[slices_1patch(i,0).ID()].push_back(i);
  }

  // sync collectively vars on each sphere
  for(int i = 0; i<nslices; ++i) {
    // sync expects a pointer to an array, standard guarantees that vector is stored contiguously in memory
    CCTK_INT* vars_ptr = &vars[i][0];
    // sync checks that all vars on the same sphere
    SphericalIntegrator_CollectiveSync(cctkGH,vars_ptr,vars[i].size());
  }

  // start to integrate each var (on timelvl 0)
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    if(cctk_iteration % slices_1patch(i,0).integrate_every() == 0) {
      // try to get output scalar index, check if it is actually there
      CCTK_INT output_index = CCTK_VarIndex(slices_1patch(i,0).outname().c_str());
      if(output_index < 0)
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "couldn't get index of output variable '%s'", slices_1patch(i, 0).outname().c_str());
      // get the pointer to the output scalar
      CCTK_REAL* result = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,output_index,NULL);
      // integrate and store result
      *result = slices_1patch(i, 0).integrate();
    }
  }
}
