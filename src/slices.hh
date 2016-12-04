

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

#ifndef _SLICES_
#define _SLICES_

#include <cctk.h>
#include <vector>
#include <utils.hh>
#include "mpi.h"
#include "spheredata_1patch.hh"

namespace SPI {


using namespace std;


/// a new radius for the slices that don't exactly lie on Llama radial-gridpoints
/// but not insist of sticking to the given radius so that we can shift the sphere radius
/// to the closest available Llama radial-gridpoint.
extern vector<CCTK_REAL> radius_internal;

/// new angular resolution in case we have Llama activated so that we can directly
/// take integer multiples of the Llama angular resolution
extern vector<int> ntheta_internal;
extern vector<int> nphi_internal;

/// a vector that stores for each slice-no the pointer to the radius storage.
extern vector<void*> radius_pointers;


/**
   This carries all data of all slices of a given type
   and assignes groups of processors that can be used for
   the various slices to get a good load balance.
*/
template <class SD>
class slices
{
   public :
            slices() : _slice(0), proc_load(0)
            { }
/*
            slices(const vector<SD>& slice_) : _slice(slice_), proc_load(0)
            { }
*/

            virtual ~slices() { }
            /// access "i-th" slices at timelevel "tl" stored in this class.
            /// This function is used to access slices (and its functions) from other thorns (if needed)
            SD operator()(const int i, const int tl) const { return _slice[i][tl]; }

            /// modify "i-th" slices at timelevel "tl" stored in this class
            /// This function is used to modify slices from other thorns (if needed)
            SD& operator()(const int i, const int tl) { return _slice[i][tl]; }

	          /// access to registered slices
            const vector<deque<SD > >& slice() const { return _slice; }

            /// convert a slice to a standard spherical surface with 1 patch by using
            /// the given number of gridpoints
            void convert_to_1patch(int const ntheta, int const nphi, CCTK_REAL* const array) const { }

            /// create storage for a new slice with some timelevels as decribed by the n-th parameter in the parfile
            /// and return the slice-id
            int register_slice(const string& varname,
                               const string& result,
                               int const slice_parameter_no,
                               int const integrate_every,
                               int const interpolate_every,
                               const integration_t integration_type,
                               const distrib_method_t distrib_method,
                               const int internal_gf_index);

            /// shifts all timelevels of i-th slice backwards, deletes the last one and creates storage for the first one
            void cycle_timelevels(const int i)
            {
               // only cycle if we have more than one timelevel
               if (_slice[i].size() > 1)
               {
                  // create new slice by using parameters of the most recent timelevel
                  _slice[i].push_front(SD(_slice[i].front().varname(), _slice[i].front().outname(), _slice[i].front().ID(),
                                          _slice[i].front().npoints()[0], _slice[i].front().npoints()[1], _slice[i].front().nghosts(),
                                          _slice[i].front().radius(0, 0, 0),
                                          _slice[i].front().radius_pointer(),
                                          _slice[i].front().origin(),
                                          _slice[i].front().has_constant_radius(),
                                          _slice[i].front().symmetry(),
                                          _slice[i].front().integrate_every(),
                                          _slice[i].front().interpolate_every(),
                                          _slice[i].front().integration_type(),
                                          _slice[i].front().distrib_method(),
                                          _slice[i].front().internal_gf_index(),
                                          _slice[i].front().processors()));
                  // remove last timelevel
                  _slice[i].pop_back();
               }
            }
   private :
            /// get a group of processors that can be used for distribution
            /// This is really only important if we intend to distribute some
            /// slice on only one processor (distrib_method == single). 
            /// In that case we want to make sure that the slices are carried
            /// by different processors
            vector<int> get_processors(const int sn, const distrib_method_t distrib_method)
            {
               DECLARE_CCTK_PARAMETERS
               // return, if no procs are needed (i.e. vol integral)
               if(distrib_method == undefined_distrib)
                return vector<int>();

               int nprocs = 1;
               MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );

               if (proc_load.size() != nprocs)
                  proc_load.resize(nprocs, 0);

               if (distrib_method == single)
               {
                  // find processor with lowest load
                  int proc = 0;
                  for (int i=0; i < nprocs; ++i)
                     if (proc_load[i] < proc_load[proc])
                        proc = i;

                  // we are going to assign the slice to processor 'proc'
                  // so we add the number of gridpoints to the load.
                  proc_load[proc] += SD::npatches*ntheta_internal[sn]*nphi_internal[sn];

                  return vector<int>(1, proc);
               }
               else
               {
                  // allow to use all available processors
                  vector<int> procs = vector<int>(nprocs, 0);
                  for (int i=0; i < nprocs; ++i)
                     procs[i] = i;

                  return procs;
               }
            }

            /// the data of all sliced variables
            vector<  // per sliced variable
                   deque<SD> > _slice;  // per timelevel

            /// stores how many spherical gridpoints each of the processors
            /// carries.
            vector<int> proc_load;
};





/// all slices for 1patch systems.
extern slices<spheredata_1patch<CCTK_REAL> > slices_1patch;

extern slices<spheredata_1patch<CCTK_REAL> > radius_1patch;

/// some shortcuts
typedef spheredata_1patch<CCTK_REAL>::integrator integrator_1patch;
typedef spheredata_1patch<CCTK_REAL>::const_iter const_iter_1patch;
typedef spheredata_1patch<CCTK_REAL>::iter       iter_1patch;

#define INDEX1P(x) x

/// given a variable-number we check if this is a 1-patch slice
inline bool is_1patch(CCTK_INT varno)
{
   if (varno >= 0 && varno < 10000)
      return true;
   return false;
}


}


#endif

