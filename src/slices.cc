
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

#include <cctk_Parameters.h>

#include "slices.hh"

namespace SPI {

template <class SD>
int slices<SD>::register_slice(const string& varname,
                               const string& outname,
                               const int slice_parameter_no,
                               const int integrate_every,
                               const int interpolate_every,
                               const integration_t integration_type,
                               const distrib_method_t distrib_method,
                               const int internal_gf_index) {
               DECLARE_CCTK_PARAMETERS
               // shortcut
               const int sn = slice_parameter_no;

               assert(slice_parameter_no < nslices);

               _slice.resize(_slice.size()+1);

               // number of current slice
               const int cs = _slice.size()-1;

               vector<int> processors = get_processors(sn, distrib_method);

               void* radii = NULL; // pointer to the data of radius function defined for the slice.
               if (varname != "ss_radius")
               {
                  // get pointer to slice that contains associated radius-function
		              assert(sn < radius_pointers.size());
                  radii = radius_pointers[sn];
               }

               // set up new slice, get storage and distribute it
               // over a defined group of processors
               SD sd(varname, outname, sn,
                     ntheta_internal[sn], nphi_internal[sn], nghostzones[sn],
                     radius_internal[sn], radii,
                     vect<CCTK_REAL,3>(origin_x[sn], origin_y[sn], origin_z[sn]),
                     set_spherical[sn],
                     vect<bool,3>(false),  // don't consider symmetries for now.
                     integrate_every,
                     interpolate_every,
                     integration_type,
                     distrib_method,
                     internal_gf_index,
                     processors);

               _slice[cs].push_back(sd);

               return _slice.size()-1;
            }

  // Instantiate template
  template class slices<spheredata_1patch<CCTK_REAL> >;

} // namespace SPI
