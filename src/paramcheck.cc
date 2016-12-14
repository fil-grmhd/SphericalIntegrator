
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



extern "C" void SphericalIntegrator_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  for(int i = 0; i<nslices; ++i) {
    if(!bnstracker_positions) {
      // check wrong activation of tracking
      if((track_star1[i]) || (track_star2[i]) ||
         (track_star1_md[i]) || (track_star2_md[i]))
        CCTK_WARN(0,"Please activate bnstracker_positions to be able to track the stars.");
    }
    else {
      if(((track_star1[i] ? 1:0) + (track_star2[i] ? 1:0)
         +(track_star1_md[i] ? 1:0) + (track_star2_md[i] ? 1:0)) > 1)
         CCTK_WARN(0,"Can't track two different positions with the same sphere!");
    }
  }
}
