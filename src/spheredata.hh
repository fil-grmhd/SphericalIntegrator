
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

#ifndef _SPHEREDATA_
#define _SPHEREDATA_

#include <cctk.h>
#include <cctk_Functions.h>
#include <util_Table.h>

#include <cstdlib>
#include <vector>


/// use Carpet-vectors for fast and easy vector-handling
#include "vect.hh"
#include "IO/io.h"
#include "int_coeffs.hh"
#include "sYlm/sYlm.h"
#include "mpi.h"


namespace SPI {


/// each slice will be initialized with this value
#define POISON_VAL 10e8


// interpolation parameters
#define N_DIMS 3
#define N_INPUT_ARRAYS 1
#define N_OUTPUT_ARRAYS 1


/// interpolation order to be used
extern int interpolator_order;


static const CCTK_REAL PI = 4.0*atan(1.0);


// distribution method: constant:           all processors hold all the data
//                      split:              the data is split across multiple processors
//                      single:             only one processor holds all the data
//                      undefined_distrib:  not defined
enum distrib_method_t { constant, split, single, undefined_distrib };

// integration type:    volume:              volume integral in sphere
//                      surface:             surface integral on sphere
//                      undefined_integral:  not defined
enum integration_t { volume, surface, undefined_integral };


/// conversion to C++ vector
template<typename T, int D>
inline vector<T>& operator<<(vector<T>& a, const vect<T, D>& b)
{
   a.resize(D);
   for (int d=0; d < D; d++)
      a[d] = b[d];

   return a;
}


// translate an absolute index to a multi-index
template <int D>
inline vect<int, D> multi_index(const int i, const vect<int, D>& num_gp)
{
   vect<int, D> multi;

   int mul = 1;
   for (int j=0; j < D-1; j++)
   {
      multi[j] = (i/mul) % num_gp[j];
      mul *= num_gp[j];
   }
   multi[D-1] = i/mul;

   return multi;
}


/**
   Abstract base class that defines the interface for a container
   that carries spherical surface data.
*/
template <typename T>
class spheredata
{
   public :
            static int const npatches = 0;
            class iter;  // forward declaration
            class const_iter;  // forward declaration

            spheredata() { }
            spheredata(const string& varname_,
                       const string& outname_,
                       const int id_,
                       const int ntheta_,
                       const int nphi_,
                       const int nghosts_,
                       const CCTK_REAL radius_,
                       void* radii_,
                       const vect<CCTK_REAL, 3>& origin_,
                       const bool has_constant_radius_,
                       const vect<bool, 3>& symmetry_,
                       const int integrate_every_,
                       const int interpolate_every_,
                       const integration_t integration_type_,
                       const distrib_method_t distrib_method_,
                       const int internal_gf_index_,
                       const vector<int>& processors_)
               : _varname(varname_), _outname(outname_), _id(id_), _name(varname_),
                 _ntheta(ntheta_+2*nghosts_), _nphi(nphi_+2*nghosts_), _nghosts(nghosts_),
                 _radius(radius_), _origin(origin_),
                 _has_constant_radius(has_constant_radius_),
                 _symmetry(symmetry_),
                 _integrate_every(integrate_every_),
                 _interpolate_every(interpolate_every_),
                 _integration_type(integration_type_),
                 _distrib_method(distrib_method_),
                 _internal_gf_index(internal_gf_index_),
                 _processors(processors_)
            {
               stringstream str;
               str << "slice=" << id_ << "::" << _name;
               _name = string(str.str());
            }

            spheredata(const spheredata& sd) { *this = sd; }
            virtual ~spheredata() { }

            /// access local surface data on patch "p", index i,j
            T operator()(const int p, const int i, const int j) const { return T(0); }

            /// modify local surface data on patch "p", index i,j
            T& operator()(const int p, const int i, const int j) { }

            /// return pointer to surface data
            void* data_pointer() const { return NULL; }

            /// access local surface radius on patch "p", index i,j
            CCTK_REAL radius(const int p, const int i, const int j) const { return T(0); }

            /// modify local surface radius on patch "p", index i,j
            CCTK_REAL& radius(const int p, const int i, const int j) { }

            /// return pointer to surface radius data
            void* radius_pointer() const { return NULL; }

            /// returns x-coordinate value of local point i,j on patch p
            CCTK_REAL cart_x(const int p, const int i, const int j) const { return 0; }
            /// returns y-coordinate value of local point i,j on patch p
            CCTK_REAL cart_y(const int p, const int i, const int j) const { return 0; }
            /// returns z-coordinate value of local point i,j on patch p
            CCTK_REAL cart_z(const int p, const int i, const int j) const { return 0; }

            /// access delta-spacing
            vect<CCTK_REAL, 2> delta() const { return 0; }

            /// access local angular coordinates (e.g. six-patch coordinate system)
            vect<CCTK_REAL, 2> coord(const int p, const int i, const int j) const { return 0; }

            /// access global angular coordinates (standard theta,phi spherical coordinate system)
            vect<CCTK_REAL, 2> coord_spherical(const int p, const int i, const int j) const { return 0; }

            /// access origin data
            vect<CCTK_REAL, 3> origin() const { return _origin; }

            /// modify origin data
            vect<CCTK_REAL, 3>& origin() { return _origin; }

            /// global surface size on patch "p" (== npoints == vect<int, 2>(ntheta, nphi))
            vect<int, 2> gsh(const int p) const { return 0; }

            /// local size on patch "p"
            vect<int, 2> lsh(const int p) const { return 0; }

            /// upper local bound on patch "p"
            vect<int, 2> ubnd(const int p) const { return 0; }

            /// upper local bound on patch "p"
            vect<int, 2> lbnd(const int p) const { return 0; }

            /// the number of global gridpoints on one patch (is supposed to be the same on all patches)
            vect<int, 2> npoints() const { return vect<int,2>(_ntheta, _nphi); }

            /// returns the global number of points (ghostpoints inclusive)
            /// on one patch (which is always the same for all patches)
            int ntheta() const { return _ntheta; }
            int nphi() const { return _nphi; }

            /// number of ghostpoints (interpatch and interprocess for all directions)
            int nghosts() const { return _nghosts; }

            /// given two surface indeices this will return the linear index
            int SINDEX2D(const int p, const int i, const int j) const
            {
               return 0;
            }

            /// returns the MPI-proc-id of the process that owns the local data
            int proc_id() const { return _proc_id; }

            /// access any theta/phi coordinate via interpolation
            CCTK_REAL interpolate(const CCTK_REAL theta, const CCTK_REAL phi) const { return 0; }

            /// interpolate from Cactus gridfunctions onto sphere
            void interpolate(const cGH* const cctkGH) { }

            /// surface integral over surface with optional function pointer
            /// that is supposed to be multiplied to each value on the sphere
            CCTK_REAL integrate(/*const CCTK_REAL* (f)(const CCTK_REAL theta, const CCTK_REAL phi)*/) const { return 0; }

            /// take pointwise derivative on patch p and point i,j in theta direction
            CCTK_REAL dx(const int p, const int i, const int j) const
            {
                return 0;
            }

            /// take pointwise derivative on patch p and point i,j in phi direction
            CCTK_REAL dy(const int p, const int i, const int j) const
            {
                return 0;
            }

            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dxdx(const int p, const int i, const int j) const
            {
                return 0;
            }

            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dxdy(const int p, const int i, const int j) const 
            {
                return 0;
            }

            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dydy(const int p, const int i, const int j) const 
            {
                return 0;
            }

            /// returns the L2-norm over the sphere
            CCTK_REAL normL2() const
            {
               return 0;
            }

            /// returns the infinity-norm over the sphere
            CCTK_REAL normLinf() const
            {
               return 0;
            }

            /// project onto spin-weighted spherical harmonic with spin
            /// s, and usual l, m
            CCTK_COMPLEX contract_with_sYlm(const int s, const int l, const int m) const { return CCTK_Cmplx(0,0); };

            /// returns the surface determinant at current point
            CCTK_REAL det(const_iter& i) const { return 0; }

            /// this flag is used to distinguish between constant spheres that
            /// don't need to store a pointwise radius
            /// and spheres that need to.
            bool has_constant_radius() const { return _has_constant_radius; }

            /// the symmetry of the sphere (symmetric_x, symmetric_y, symmetric_z)
            vect<bool, 3> symmetry() const { return _symmetry; }

            /// returns the name of the slice
            string name() const { return _name; }

            /// returns the original name of the Cactus variable that is sliced
            string varname() const { return _varname; }

            /// returns the original name of the Cactus scalar that gets the integration result
            string outname() const { return _outname; }

            /// returns pointer to output variable, defined by outname()
            CCTK_REAL* outpointer(const cGH* const cctkGH) const {
              // try to get output scalar index, check if it is actually there
              CCTK_INT output_index = CCTK_VarIndex(this->outname().c_str());
              if(output_index < 0)
                CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "Couldn't get index of output variable '%s' for '%s' on sphere %i",
                           this->outname().c_str(),
                           this->varname().c_str(),
                           this->ID());
              // get the pointer to the output scalar
              return (CCTK_REAL*) CCTK_VarDataPtrB(cctkGH,0,output_index,NULL);

            }

            /// returns integration iteration number
            int integrate_every() const { return _integrate_every; }

            /// returns interpolation iteration number
            int interpolate_every() const { return _interpolate_every; }

            /// returns internal gf index, for volume integrals
            int internal_gf_index() const { return _internal_gf_index; }

            /// returns the ID of the slice == n-th spherical slice in parfile
            int ID() const { return _id; }

            /// return decomposed local piece of data for this MPI process
            void decompose() { }

            /// returns the distribution method
            distrib_method_t distrib_method() const { return _distrib_method; }

            /// returns the integration type
            integration_t integration_type() const { return _integration_type; }

            /// a group of processors that carry parts of the slice
            vector<int> processors() const { return _processors; }

            /// write the SphereData including all attributes
            /// to the output stream defined by io_base (or its inheritants)
            void output (io_base& io) const { }

            void input (io_base& io) { }

            iter begin() { return iter(); }

            /// iterator class to traverse the grid of the slice
            class const_iter
            {
               public :
                        const_iter() { }
                        virtual ~const_iter() { }

                        /// dereferencing operator as data-point accessor
                        T operator*() const { return 0; }

                        spheredata<T> const * to_spheredata() const { return NULL; }

                        void operator++() { (*this)++; }

                        void operator++(int) { }

                        /// query whether iterator is done
                        bool done() const
                        {
                           return false;
                        }

                        /// index-struct
                        struct idx_t
                        {
                           /// current patch
                           int p;
                           /// current grid indices
                           int i, j;
                           /// curent lienar index
                           int ij;
                           /// current local coordinates
                           CCTK_REAL theta, phi;
                        };

                        idx_t idx() const { return _idx; }

               protected :
                        idx_t _idx;
            };

            class iter : public const_iter 
            {
               public :
                        iter() : const_iter() { }
                        virtual ~iter() { }

                        /// dereferencing operator as data-point accessor
                        // T& operator*() { return (*this->_spheredata).data[this->_idx.ij];  }
                        T& operator*() { abort(); }

                        spheredata<T>* to_spheredata()
                        {
                           return NULL; 
                        }
            };

            class integrator
            {
               public :
                        integrator() { }
                        virtual ~integrator() { }

                        void init() { }
                        void finalize() { }
                        void sum(const_iter& it, CCTK_REAL f = 1.0, CCTK_REAL det = 1.0) { }
            };


   protected :

            /// for now we will not work with symmetries....
            vect<bool, 3> _symmetry;

            bool _has_constant_radius;

            /// the distribution method
            distrib_method_t _distrib_method;

            /// the integration type
            integration_t _integration_type;

            ///...and the corresponding process that carries the data
            int _proc_id;

            /// name of the slice (sliced variable-name)
            string _name;

            /// the id of this slice
            int _id;

            /// the name of the Cactus gridfunction which we want to slice
            string _varname;

            /// the name of the Cactus scalar, where the integration result should be stored
            string _outname;

            /// how often should the GF be integrated?
            int _integrate_every;

            /// how often should the GF be interpolated?
            int _interpolate_every;

            /// index of the internal gf, used for volume integration
            int _internal_gf_index;

            /// a constant radius
            CCTK_REAL _radius;

            /// the slice's origin
            vect<CCTK_REAL, 3> _origin;

            /// the resolution of one patch
            int _ntheta, _nphi;

            /// the number of ghostpoints for each patch and direction
            int _nghosts;

            /// the processor-ids among which to distribute the slice
            vector<int> _processors;

            /// a flag specifying whether this surface is valid or not.
            bool _valid;
};



template <typename T>
inline io_base& operator<< (io_base& io, const spheredata<T>& sd)
{
   sd.output(io);

   return io;
}


}


#endif

