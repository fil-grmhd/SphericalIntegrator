
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
#include "carpet.hh"

#include <cmath>

using namespace SPI;



extern "C" CCTK_INT SphericalIntegrator_SurfaceSync(const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_INT varno)
{
   DECLARE_CCTK_PARAMETERS

   // CHECK WTF?
   assert(Carpet::do_global_mode);

   cGH const * restrict const cctkGH = (cGH*) cctkGH_;

   assert(varno >= 0);

   assert(varno < slices_1patch.slice().size());
   slices_1patch.cycle_timelevels(varno);
   slices_1patch(varno, 0).interpolate(cctkGH);

   return 0;
}

extern "C" CCTK_INT SphericalIntegrator_CollectiveSurfaceSync(const CCTK_POINTER_TO_CONST cctkGH_, const CCTK_POINTER_TO_CONST varno_, const CCTK_INT number_of_vars)
{
   DECLARE_CCTK_PARAMETERS

   assert(Carpet::do_global_mode);

   cGH const * restrict const cctkGH = (cGH*) cctkGH_;
   int const * const varno = (int*) varno_;

   if (number_of_vars == 0) return 0;

   assert(varno != NULL);

   if (!use_carpet_interp1)
   {
      // make sure all variable numbers belong to the same slice!
      const int ID = slices_1patch(varno[0], 0).ID();
      for (int i=1; i < number_of_vars; ++i)
      {
         assert(ID == slices_1patch(varno[i], 0).ID());
         assert(varno[i] < slices_1patch.slice().size());
         slices_1patch.cycle_timelevels(varno[i]);
      }

      // get the input variable's index
      vector<int> varindices (number_of_vars, -1);

      for (int i=0; i < number_of_vars; ++i)
      {
         varindices[i] = CCTK_VarIndex(slices_1patch(varno[i], 0).varname().c_str());
         if (varindices[i] < 0)
            CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "couldn't get index of slice variable '%s'", slices_1patch(varno[i], 0).varname().c_str());
      }

      vector<CCTK_REAL*> values (number_of_vars, static_cast<CCTK_REAL*>(NULL));

      for (int i=0; i < number_of_vars; ++i)
         values[i] = (CCTK_REAL*) slices_1patch(varno[i], 0).data_pointer();

      // set up the interpolation
      assert(interp_setups.size() > slices_1patch(varno[0], 0).ID());
      interp_setup_t* &interp_setup = interp_setups[slices_1patch(varno[0], 0).ID()];
      if (not interp_setup or interp_setup->fasterp_setup->outofdate()) {
         if (interp_setup)
           delete interp_setup;
         interp_setup = new interp_setup_t(slices_1patch(varno[0], 0).lsh(0)[0] * slices_1patch(varno[0], 0).lsh(0)[1]);

         // allocate storage for coordinates
         fasterp_glocs_t locations (interp_setup->npoints);

         // get Cartesian coordinate values of gridpoints on spherical surface
         for (int j=0; j < slices_1patch(varno[0], 0).lsh(0)[1]; ++j) {
         for (int k=0; k < slices_1patch(varno[0], 0).lsh(0)[0]; ++k) {
            const int l = k + slices_1patch(varno[0], 0).lsh(0)[0]*j;
            double x = slices_1patch(varno[0], 0).cart_x(0, k, j);
            double y = slices_1patch(varno[0], 0).cart_y(0, k, j);
            double z = slices_1patch(varno[0], 0).cart_z(0, k, j);
            /*if (have_symmetries)
            {
               if (x < 0)
               {
                  x = -x;
                  y = -y;
               }
               if (y < 0)
               {
                  double const t = x;
                  x = -y;
                  y =  t;
               }
               if (z < 0)
               {
                  z = -z;
               }
            }*/
            locations.coords[0][l] = x;
            locations.coords[1][l] = y;
            locations.coords[2][l] = z;
         }
         }

         interp_setup->fasterp_setup =
         new fasterp_setup_t(cctkGH, locations, interpolator_order);
      }

      // do the interpolation
      assert(interp_setup->fasterp_setup);
      assert(interp_setup->npoints == slices_1patch(varno[0], 0).lsh(0)[0] * slices_1patch(varno[0], 0).lsh(0)[1]);
      interp_setup->fasterp_setup->interpolate (cctkGH, varindices, values);

      if (not slices_1patch(varno[0], 0).has_constant_radius()) {
         delete interp_setup;
         interp_setup = NULL;
      }
   }

   if(use_carpet_interp1)
   {
      // make sure all variable numbers belong to the same slice!
      const int ID = slices_1patch(varno[0], 0).ID();
      for (int i=1; i < number_of_vars; ++i)
      {
         assert(ID == slices_1patch(varno[i], 0).ID());
         assert(varno[i] < slices_1patch.slice().size());
         slices_1patch.cycle_timelevels(varno[i]);
      }

      const void*   interp_coords[N_DIMS];
      vector<CCTK_INT> input_array_indices(number_of_vars);
      vector<void*> output_arrays(number_of_vars);
      const vector<CCTK_INT> input_array_type_codes(number_of_vars, CCTK_VARIABLE_REAL);
               //= { CCTK_VARIABLE_REAL };
      const vector<CCTK_INT> output_array_type_codes(number_of_vars, CCTK_VARIABLE_REAL);
               //= { CCTK_VARIABLE_REAL };

      // set interpolation handles
      int operator_handle = CCTK_InterpHandle("Lagrange polynomial interpolation");
      if (operator_handle < 0)
         CCTK_WARN(0, "can’t get interpolation handle!");
      int param_table_handle = Util_TableCreateFromString("order=4");
      if (param_table_handle < 0)
         CCTK_WARN(0, "can’t create parameter table!");
      int coordsys_handle = CCTK_CoordSystemHandle("cart3d");
      if (coordsys_handle < 0)
         CCTK_WARN(0, "can’t create coordsys handle! Forgot to activate a Coordinate-Thorn?");

      for (int i=0; i < number_of_vars; ++i)
      {
         input_array_indices[i] = CCTK_VarIndex(slices_1patch(varno[i], 0).varname().c_str());
         if (input_array_indices[i] < 0)
            CCTK_WARN(0, "error getting VarIndex of variable that shall be sliced");
      }

     // project variable onto sphere
     // this means we have to 3d-interpolate every non-coinciding gridpoint of
     // the spherical grid from neighboring gridpoints of the Cartesian grid.
     // for this, we first calculate the Cartesian coordinate of the jk-th spherical
     // surface gridpoint
     int n_interp_points = slices_1patch(varno[0], 0).lsh(0)[0]
                          *slices_1patch(varno[0], 0).lsh(0)[1];
                          //lsh(m)[0]*lsh(m)[1];

     // Arrays of Cartesian coordinates of the surface points onto which we want to interpolate
     vector<CCTK_REAL> interp_x(n_interp_points, POISON_VAL);
     vector<CCTK_REAL> interp_y(n_interp_points, POISON_VAL);
     vector<CCTK_REAL> interp_z(n_interp_points, POISON_VAL);

     // get Cartesian coordinate values of gridpoints on spherical surface
     for (int j=0; j < slices_1patch(varno[0], 0).lsh(0)[1]; ++j)
        for (int k=0; k < slices_1patch(varno[0], 0).lsh(0)[0]; ++k)
        {
           const int l = k + slices_1patch(varno[0], 0).lsh(0)[0]*j;
           interp_x[l] = slices_1patch(varno[0], 0).cart_x(0, k, j);
           interp_y[l] = slices_1patch(varno[0], 0).cart_y(0, k, j);
           interp_z[l] = slices_1patch(varno[0], 0).cart_z(0, k, j);
        }

     interp_coords[0] = (const void*) &interp_x.front();
     interp_coords[1] = (const void*) &interp_y.front();
     interp_coords[2] = (const void*) &interp_z.front();

     for (int i=0; i < number_of_vars; ++i)
        output_arrays[i] = slices_1patch(varno[i], 0).data_pointer();

     // Do the actual interpolation.
     // Only those processes interpolate that contain data
     // All other processes simply send their
     // grid-data.
     //if (data[m].size() == 0)
     //   n_interp_points = 0;      // all other processes shall not interpolate

     if (CCTK_InterpGridArrays(cctkGH,
                                N_DIMS,
                                operator_handle, param_table_handle,
                                coordsys_handle,
                                n_interp_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                number_of_vars,
                                   &input_array_indices.front(),
                                number_of_vars,
                                   &output_array_type_codes.front(),
                                   &output_arrays.front()) < 0)
        CCTK_WARN(1, "error return from interpolator!");
   }

   return 0;
}


extern "C" void SphericalIntegrator_CollectiveInterpolation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  vector<CCTK_INT> vars[nslices];

  if(verbose > 0)
    CCTK_INFO("Interpolating surface variables to spheres.");

  // collect all vars and sort them to their sphere ID
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    // only interpolate if it is wanted
    if((slices_1patch(i,0).interpolate_every() != 0)
      && (cctk_iteration % slices_1patch(i,0).interpolate_every() == 0))
        vars[slices_1patch(i,0).ID()].push_back(i);
  }

  // sync collectively vars on each sphere
  for(int i = 0; i<nslices; ++i) {
    // sync expects a pointer to an array
    // (standard guarantees that vector is stored contiguously in memory)
    CCTK_INT* vars_ptr = &vars[i][0];
    // sync checks that all vars on the same sphere
    SphericalIntegrator_CollectiveSurfaceSync(cctkGH,vars_ptr,vars[i].size());
  }
}

extern "C" void SphericalIntegrator_CollectiveVolumeSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose > 0)
    CCTK_INFO("Syncing volume variables to internal GFs and mask by spheres.");

  // vector of volume integrated vars
  vector<CCTK_INT> vol_vars;

  // collect volume integration vars
  for(int i = 0; i<slices_1patch.slice().size(); ++i) {
    if((slices_1patch(i,0).integration_type() == volume)
        && (slices_1patch(i,0).integrate_every() != 0)
        && (cctk_iteration % slices_1patch(i,0).integrate_every() == 0)) {
      if(slices_1patch(i,0).has_constant_radius())
        vol_vars.push_back(i);
      else
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Volume integrals over elliptical surfaces not supported yet. (sn = %i)",slices_1patch(i,0).ID());
    }
  }
  // drop if nothing is integrated this iteration
  if(vol_vars.size() == 0)
    return;

  // check if enough tmp GFs are available
  if(vol_vars.size() > max_volume_integrals)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Please increase max_volume_integrals parameter to at least %lu", vol_vars.size());

  // get index size of GFs
  const CCTK_INT gf_size = UTILS_GFSIZE(cctkGH);

  // vector of pointers to the actual and temporary GFs
  vector<CCTK_REAL*> vol_vars_pointers(vol_vars.size(),NULL);
  vector<CCTK_REAL*> vol_vars_tmp_pointers(vol_vars.size(),NULL);

  for(int i = 0; i<vol_vars.size(); ++i) {
    // get pointer to internal tmp gf
    CCTK_REAL* tmp_gf_pointer = &ss_tmp_volume_gfs[slices_1patch(vol_vars[i],0).internal_gf_index()*gf_size];
    // store in vector for loop below
    vol_vars_tmp_pointers[i] = tmp_gf_pointer;

    // get pointer to actual GF
    CCTK_INT gf_index = CCTK_VarIndex(slices_1patch(vol_vars[i],0).varname().c_str());
    if(gf_index < 0)
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "couldn't get index of variable '%s'", slices_1patch(vol_vars[i], 0).varname().c_str());
    CCTK_REAL* gf_pointer = (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,gf_index,NULL);
    // store it in vector for loop below
    vol_vars_pointers[i] = gf_pointer;
  }

  #pragma omp parallel for schedule(static)
  for(int ijk = 0; ijk < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; ++ijk) {
    // compute det(metric) for volume element
    CCTK_REAL det_g = utils::metric::spatial_det(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk]);
    for(int i = 0; i<vol_vars.size(); ++i) {
      // calculate (coordinate) distance from sphere origin
      CCTK_REAL x_dist = x[ijk]-slices_1patch(vol_vars[i],0).origin()[0];
      CCTK_REAL y_dist = y[ijk]-slices_1patch(vol_vars[i],0).origin()[1];
      CCTK_REAL z_dist = z[ijk]-slices_1patch(vol_vars[i],0).origin()[2];
      CCTK_REAL distance = std::sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist);

      // set to zero outside of the sphere
      if(distance <= slices_1patch(vol_vars[i],0).radius())
        vol_vars_tmp_pointers[i][ijk] = det_g*vol_vars_pointers[i][ijk];
      else
        vol_vars_tmp_pointers[i][ijk] = 0.0;
    }
  }
}
