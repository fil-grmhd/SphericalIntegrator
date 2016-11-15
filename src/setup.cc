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

#include "setup.hh"



namespace SPI {

static
CCTK_REAL min (CCTK_REAL const x, CCTK_REAL const y)
{
  return x < y ? x : y;
}

static
CCTK_REAL max (CCTK_REAL const x, CCTK_REAL const y)
{
  return x > y ? x : y;
}

#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

slices<spheredata_1patch<CCTK_REAL> > slices_1patch;

slices<spheredata_1patch<CCTK_REAL> > radius_1patch;


/// a flag for each of the slices defining whether it can take advantage of Llama
vector<bool> can_use_Llama_internal;

/// a new radius for the slices that don't exactly lie on Llama radial-gridpoints
/// but not insist of sticking to the given radius so that we can shift the sphere radius
/// to the closest available Llama radial-gridpoint.
vector<CCTK_REAL> radius_internal;

/// new angular resolution in case we have Llama activated so that we can directly
/// take integer multiples of the Llama angular resolution
vector<int> ntheta_internal;
vector<int> nphi_internal;

/// a vector that stores for each slice-no the pointer to the radius storage.
vector<void*> radius_pointers;


/// if a Multipatch system is present that supports Thornburg04-coordinates,
/// the Llama gets activated
bool Llama_activated = false;

/// interpolation order to be used
int interpolator_order;

/// interpolation setups for slices
vector<interp_setup_t *> interp_setups;
}


using namespace SPI;


extern "C" void SphericalIntegrator_Setup(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   CCTK_INFO("Basic setup of slices.");

   Llama_activated = false;

   can_use_Llama_internal = vector<bool>(nslices, false);
   radius_internal = vector<CCTK_REAL>(nslices, 0);
   ntheta_internal = vector<int>(nslices, 0);
   nphi_internal = vector<int>(nslices, 0);


   // now check each slice
   for (int i=0; i < nslices; ++i)
   {
      can_use_Llama[i] = false;
      // set new radius and resolution
      // just use given settings
      new_radius[i] = radius[i];
      ss_ntheta[i] = ntheta[i];
      ss_nphi[i] = nphi[i];

      // set ss_info
      ss_origin_x[i] = origin_x[i];
      ss_origin_y[i] = origin_y[i];
      ss_origin_z[i] = origin_z[i];

      // mark surface as unitialized
      ss_active[i] = 0;
      ss_valid[i] = 0;

      if (set_spherical[i])
      {
         ss_area[i] = 4 * pi * pow (radius[i], 2);

         ss_mean_radius[i] = new_radius[i];

         ss_centroid_x[i] = origin_x[i];
         ss_centroid_y[i] = origin_y[i];
         ss_centroid_z[i] = origin_z[i];

         ss_quadrupole_xx[i] = 0.0;
         ss_quadrupole_xy[i] = 0.0;
         ss_quadrupole_xz[i] = 0.0;
         ss_quadrupole_yy[i] = 0.0;
         ss_quadrupole_yz[i] = 0.0;
         ss_quadrupole_zz[i] = 0.0;

         ss_min_radius[i] = new_radius[i];
         ss_max_radius[i] = new_radius[i];

         ss_min_x[i] = origin_x[i] - new_radius[i];
         ss_min_y[i] = origin_y[i] - new_radius[i];
         ss_min_z[i] = origin_z[i] - new_radius[i];
         ss_max_x[i] = origin_x[i] + new_radius[i];
         ss_max_y[i] = origin_y[i] + new_radius[i];
         ss_max_z[i] = origin_z[i] + new_radius[i];

         // radii/origins are valid
         ss_active[i] = 1;
         ss_valid[i] = 1;
      }
      else if (set_elliptic[i]) 
      {
         CCTK_REAL const rx2 = pow (radius_x[i], 2);
         CCTK_REAL const ry2 = pow (radius_y[i], 2);
         CCTK_REAL const rz2 = pow (radius_z[i], 2);

         ss_active[i] = +1;
         ss_valid[i] = +1;

         ss_area[i] = 0 * 4 * pi * pow (radius[i], 2);

         ss_mean_radius[i] = 0 * radius[i];

         ss_centroid_x[i] = origin_x[i];
         ss_centroid_y[i] = origin_y[i];
         ss_centroid_z[i] = origin_z[i];

         ss_quadrupole_xx[i] = 1.0 / rz2;
         ss_quadrupole_xy[i] = 0.0;
         ss_quadrupole_xz[i] = 0.0;
         ss_quadrupole_yy[i] = 1.0 / ry2;
         ss_quadrupole_yz[i] = 0.0;
         ss_quadrupole_zz[i] = 1.0 / rx2;

         ss_min_radius[i] = min (radius_x[i], min (radius_y[i], radius_z[i]));
         ss_max_radius[i] = max (radius_x[i], max (radius_y[i], radius_z[i]));

         ss_min_x[i] = origin_x[i] - radius_x[i];
         ss_min_y[i] = origin_y[i] - radius_y[i];
         ss_min_z[i] = origin_z[i] - radius_z[i];
         ss_max_x[i] = origin_x[i] + radius_x[i];
         ss_max_y[i] = origin_y[i] + radius_y[i];
         ss_max_z[i] = origin_z[i] + radius_z[i];
      }
   }
}




extern "C" void SphericalIntegrator_PostSetup(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS

   CCTK_INFO("Setup radius of slices.");

   can_use_Llama_internal = vector<bool>(nslices, false);
   radius_internal = vector<CCTK_REAL>(nslices, 0);
   ntheta_internal = vector<int>(nslices, 0);
   nphi_internal = vector<int>(nslices, 0);

   // now check each slice
   for (int i=0; i < nslices; ++i)
   {
      ntheta_internal[i] = ss_ntheta[i];
      nphi_internal[i] = ss_nphi[i];
      can_use_Llama_internal[i] = can_use_Llama[i];
      radius_internal[i] = new_radius[i];

      // We are now going to set up the radii as separately registered slice for each slice
      // Each processor should have the radius available. That means we gonna use a constant distribution.
      //ss_radius_id[i] = SphericalIntegrator_Register("ss_radius", i, 1, "const");
      ss_radius_id[i] = ONEPATCH_SLICE_IDS + radius_1patch.register_slice("ss_radius", i, 1, constant);

      if (set_elliptic[i])
      {
         CCTK_REAL const rx2 = pow (radius_x[i], 2);
         CCTK_REAL const ry2 = pow (radius_y[i], 2);
         CCTK_REAL const rz2 = pow (radius_z[i], 2);

        for (iter_1patch it=radius_1patch(INDEX1P(ss_radius_id[i]), 0).begin(); !it.done(); ++it)
        {
            CCTK_REAL const theta = it.idx().theta;
            CCTK_REAL const phi = it.idx().phi;
            CCTK_REAL const x2 = pow (sin(theta) * cos(phi), 2);
            CCTK_REAL const y2 = pow (sin(theta) * sin(phi), 2);
            CCTK_REAL const z2 = pow (cos(theta)           , 2);
            *it = 1.0 / sqrt (x2 / rx2 + y2 / ry2 + z2 / rz2);
        }
      }
   }

   // after all radius-slices have been set, we need to update
   // the pointers.
   // From this point on, radius_1patch is NOT allowed to change, otherwise
   // the radius-pointers get screwed!!!
   radius_pointers = vector<void*>(nslices, static_cast<void*>(NULL));
   for (int i=0; i < nslices; ++i)
   {
      // get radius pointer
      radius_pointers[i] = radius_1patch(INDEX1P(ss_radius_id[i]), 0).radius_pointer();
   }
   SPI::interpolator_order = interpolator_order;
   interp_setups.resize(nslices, NULL);
}
