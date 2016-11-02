
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



namespace SPS {

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
slices<spheredata_2patch<CCTK_REAL> > slices_2patch;
slices<spheredata_6patch<CCTK_REAL> > slices_6patch;

slices<spheredata_1patch<CCTK_REAL> > radius_1patch;
slices<spheredata_2patch<CCTK_REAL> > radius_2patch;
slices<spheredata_6patch<CCTK_REAL> > radius_6patch;


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


using namespace SPS;


extern "C" void SphericalSlice_Setup(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_PARAMETERS
   
   CCTK_INFO("Basic setup of slices.");
   
   // check if there is a Thornburg04 Multipatch-system available
   
   // it is not sufficient if thorn "Coordinates" is activated since
   // other thorns than Llama could also have that name.
   // Llama provides a function that just confirms that it is
   // a Thornburg04-six-patch system. 
   if (CCTK_IsFunctionAliased("MultiPatch_ProvidesThornburg04"))
   {
      if (MultiPatch_ProvidesThornburg04())
      {
         Llama_activated = true;
         CCTK_INFO("Llama sighted!");
      }
      else
         Llama_activated = false;
   }
   else
      Llama_activated = false;
   
   
   can_use_Llama_internal = vector<bool>(nslices, false);
   radius_internal = vector<CCTK_REAL>(nslices, 0);
   ntheta_internal = vector<int>(nslices, 0);
   nphi_internal = vector<int>(nslices, 0);
   
   
   // now check each slice
   for (int i=0; i < nslices; ++i)
   {
      // check, if this slice can take advantage of Llama
      if (Llama_activated)
      {
         if (CCTK_Equals(type[i], "6patch") &&
             set_spherical[i] && 
             origin_x[i] == 0 && 
             origin_y[i] == 0 && 
             origin_z[i] == 0 &&
             use_Llama[i] == 1 &&
             insist_on_radius[i] == 0 && 
             radius[i] > MultiPatch_GetInnerRadius())
         {
            can_use_Llama[i] = true;
            
         }
         else
            can_use_Llama[i] = false;
      }
      else
         can_use_Llama[i] = false;
      
      // set new radius and resolution
      if (can_use_Llama[i])
      {
         if (CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"))
         {
            vector<CCTK_REAL> physical_min(3);
            vector<CCTK_REAL> physical_max(3);
            vector<CCTK_REAL> interior_min(3);
            vector<CCTK_REAL> interior_max(3);
            vector<CCTK_REAL> exterior_min(3);
            vector<CCTK_REAL> exterior_max(3);
            vector<CCTK_REAL> spacing(3);
            
            MultiPatch_GetDomainSpecification(1, 3, &physical_min.front(), &physical_max.front(),
                                                    &interior_min.front(), &interior_max.front(),
                                                    &exterior_min.front(), &exterior_max.front(),
                                                    &spacing.front());
            
            // find closest radial gridpoint
            int r = round((radius[i] - physical_min[2]) / spacing[2]);
            new_radius[i] = r*spacing[2] + physical_min[2];
         
            ss_ntheta[i] = (physical_max[0] - physical_min[0]) / spacing[0] + 1;
            ss_nphi[i] =  (physical_max[1] - physical_min[1]) / spacing[1] + 1;
            
            if (ss_ntheta[i] % stride[i] != 0)
               CCTK_WARN(1, "The Llama is spitting on you: Not able to make use of stride since Coordinates::n_angular is not dividable by stride! Setting ntheta = Coordinates::n_angular");
            else
               ss_ntheta[i] = ss_ntheta[i] / stride[i];
            
            if (ss_nphi[i] % stride[i] != 0)
               CCTK_WARN(1, "The Llama is spitting on you: Not able to make use of stride since Coordinates::n_angular is not dividable by stride! Setting nphi = Coordinates::n_angular");
            else
               ss_nphi[i] = ss_nphi[i] / stride[i];
            
            CCTK_VInfo (CCTK_THORNSTRING, "Riding the Llama: Setting radius[%d] = %f",
                        i, new_radius[i]);
            CCTK_VInfo (CCTK_THORNSTRING, "Riding the Llama: Setting ntheta[%d] = %d",
                        i, int(ss_ntheta[i]));
            CCTK_VInfo (CCTK_THORNSTRING, "Riding the Llama: Setting nphi[%d] = %d",
                        i, int(ss_nphi[i]));
         }
         else // just use given settings
         {
            new_radius[i] = radius[i];
            ss_ntheta[i] = ntheta[i];
            ss_nphi[i] = nphi[i];
         }
         
      }
      else // just use given settings
      { 
         new_radius[i] = radius[i];
         ss_ntheta[i] = ntheta[i];
         ss_nphi[i] = nphi[i];
      }
      
      
      if (CCTK_Equals(type[i], "1patch")) 
         ss_npatches[i] = 1;
      if (CCTK_Equals(type[i], "6patch")) 
         ss_npatches[i] = 6;
      
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
         ss_active[i] = +1;
         ss_valid[i] = +1;
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




extern "C" void SphericalSlice_PostSetup(CCTK_ARGUMENTS)
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
      //ss_radius_id[i] = SphericalSlice_Register("ss_radius", i, 1, "const");
      if (CCTK_Equals(type[i], "1patch"))
         ss_radius_id[i] = ONEPATCH_SLICE_IDS + radius_1patch.register_slice("ss_radius", i, 1, constant);
      if (CCTK_Equals(type[i], "6patch"))
         ss_radius_id[i] = SIXPATCH_SLICE_IDS + radius_6patch.register_slice("ss_radius", i, 1, constant);
      
      
      
      if (set_elliptic[i]) 
      {
         CCTK_REAL const rx2 = pow (radius_x[i], 2);
         CCTK_REAL const ry2 = pow (radius_y[i], 2);
         CCTK_REAL const rz2 = pow (radius_z[i], 2);
         
         if (ss_npatches[i] == 1)
         {
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
         
         if (ss_npatches[i] == 6)
         {
            for (iter_6patch it=radius_6patch(INDEX6P(ss_radius_id[i]), 0).begin(); !it.done(); ++it)
            {
               CCTK_REAL const theta = radius_6patch(INDEX6P(ss_radius_id[i]), 0).coord_spherical(it.idx().p, it.idx().i, it.idx().j)[0];
               CCTK_REAL const phi = radius_6patch(INDEX6P(ss_radius_id[i]), 0).coord_spherical(it.idx().p, it.idx().i, it.idx().j)[1];
               CCTK_REAL const x2 = pow (sin(theta) * cos(phi), 2);
               CCTK_REAL const y2 = pow (sin(theta) * sin(phi), 2);
               CCTK_REAL const z2 = pow (cos(theta)           , 2);
               *it = 1.0 / sqrt (x2 / rx2 + y2 / ry2 + z2 / rz2);
            }
         }

      }
      
   }
   
   // after all radius-slices have been set, we need to update
   // the pointers.
   // From this point on, radius_1patch and radius_6patch are NOT allowed to change, otherwise
   // the radius-pointers get screwed!!!
   radius_pointers = vector<void*>(nslices, static_cast<void*>(NULL));
   for (int i=0; i < nslices; ++i)
   {
      // get radius pointer
      if (is_1patch(ss_radius_id[i]))
         radius_pointers[i] = radius_1patch(INDEX1P(ss_radius_id[i]), 0).radius_pointer();
      if (is_6patch(ss_radius_id[i]))
         radius_pointers[i] = radius_6patch(INDEX6P(ss_radius_id[i]), 0).radius_pointer();
   }
   
   SPS::interpolator_order = interpolator_order;
   interp_setups.resize(nslices, NULL);
}
