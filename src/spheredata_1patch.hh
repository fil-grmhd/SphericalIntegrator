
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

#ifndef _SPHEREDATA1_
#define _SPHEREDATA1_

#include <vector>

#include "mpi.h"
#include "carpetinterp2.hh"

#include "spheredata.hh"
#include "commstack.hh"

//#define DEBUG

namespace SPS {


using namespace std;
using namespace CarpetInterp2;


struct interp_setup_t {
  fasterp_setup_t * fasterp_setup;
  int npoints;

  interp_setup_t (int npoints_) : fasterp_setup (NULL), npoints (npoints_) {}
  ~interp_setup_t () { if (fasterp_setup) delete fasterp_setup; }
};
extern vector<interp_setup_t *> interp_setups;

/**
   The standard "old-school" SphericalSurface-style spherical 
   slice description.
   This will be used for simple user output (if requested)
   or other thorns that want to work on such simpler slices.
*/
template <typename T>
class spheredata_1patch : public spheredata<T>
{
   public :
            static int const npatches = 1;
            class iter;  // forward declaration
            class const_iter;  // forward declaration
            
            spheredata_1patch() { }
            spheredata_1patch(const spheredata_1patch& sd) 
            { 
               *this = sd; 
               if (spheredata<T>::varname() == "ss_radius") 
                  _radii = &data.front();
            }
            spheredata_1patch(const string& varname_,
                              const int id_,
                              const int ntheta_,
                              const int nphi_,
                              const int nghosts_,
                              const CCTK_REAL radius_,
                              void* radii_,
                              const vect<CCTK_REAL, 3>& origin_,
                              const bool has_constant_radius_,
                              const vect<bool, 3>& symmetry_,
                              const distrib_method_t distrib_method_,
                              const vector<int>& processors_,
                              const bool can_use_Llama_)
               : spheredata<T>(varname_, id_, ntheta_, nphi_, nghosts_, radius_, radii_,
                               origin_, has_constant_radius_, symmetry_,
                               distrib_method_, processors_, can_use_Llama_),
                 _lsh(vect<int,2>(ntheta_+2*nghosts_, nphi_+2*nghosts_)), _gsh(vect<int,2>(ntheta_+2*nghosts_, nphi_+2*nghosts_)),
                 _lbnd(0), _ubnd(vect<int,2>(ntheta_-1+2*nghosts_, nphi_-1+2*nghosts_)),
                 _dtheta(PI/(ntheta_)),
                 _dphi(2.0*PI/(nphi_)),
                 _radii(reinterpret_cast<CCTK_REAL*>(radii_))
            {
               _origin_theta = (-nghosts_+0.5)*_dtheta;  // theta direction is staggered
               _origin_phi   = -nghosts_*_dphi;
               
               if (spheredata<T>::_distrib_method == constant) // all processors will carry all data
               {
                  data.resize(spheredata<T>::npoints()[0]*spheredata<T>::nphi());
                  
                  /*if (!spheredata<T>::has_constant_radius())
                     _radii.resize(spheredata<T>::npoints()[0]*spheredata<T>::nphi());*/
               }
               if (spheredata<T>::_distrib_method == single) // one processor will carry all data
               {
                  int this_proc = 0;
                  assert(spheredata<T>::_processors.size() == 1);
                  MPI_Comm_rank(MPI_COMM_WORLD, &this_proc);
                  if (this_proc == spheredata<T>::_processors[0])
                  {
                     data.resize(spheredata<T>::npoints()[0]*spheredata<T>::nphi());
                  
                     /*if (!spheredata<T>::has_constant_radius())
                        _radii.resize(spheredata<T>::npoints()[0]*spheredata<T>::nphi());*/
                  }
                  else
                  {
                     _lsh = vect<int,2>(0, 0);
                  }
               }
               if (spheredata<T>::_distrib_method == split) // a group of processors will carry the data
               {
                  // do a domain decomposition
                  decompose();
                  
                  data.resize(lsh(0)[0]*lsh(0)[1]);
                  
                  /*if (!spheredata<T>::has_constant_radius())
                     _radii.resize(lsh(0)[0]*lsh(0)[1]);*/
               }
               
               // set to poison values initially
               for (int i=0; i < data.size(); ++i)
                  data[i] = POISON_VAL;
               
               // if the radius-pointer has not been set then we assume that this slice
               // cariies the radius as data!
               if (spheredata<T>::varname() == "ss_radius")//_radii == NULL)
               {
                  _radii = &data.front();
                  
                  // set radius according to initial constant value
                  for (iter i=begin(); !i.done(); ++i)
                     *i = radius_;
               }
               //if (spheredata<T>::varname() != "ss_radius")
               //   for (iter i=begin(); !i.done(); ++i)
               //      cout << _radii[i.idx().ij]<< endl;
               
               if (lsh(0)[0]-2*spheredata<T>::nghosts() > 0)
                  int_weights = vector<CCTK_REAL>(lsh(0)[0]-2*spheredata<T>::nghosts(), 0);
               
               // initialize Gauss-integration weights
               for (int i=0; i < lsh(0)[0]; ++i) 
               {
                  // dont't calculate weights in ghostzones
                  if (i < spheredata<T>::nghosts() || 
                      i > lsh(0)[0]-spheredata<T>::nghosts()-1)
                     continue;
                  
                  int_weights[i-spheredata<T>::nghosts()] = weight(i+lbnd(0)[0]-spheredata<T>::nghosts(), gsh(0)[0]-2*spheredata<T>::nghosts());
                  
               }
            }
            
            virtual ~spheredata_1patch() { }
            
            /// given two surface indices this will return the linear index
            int SINDEX2D(const int p, const int i, const int j) const
            {
               return i + lsh(p)[0]*j;
            }
            
            /// access local surface data on patch "p", index i,j
            T operator()(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
               return data[i + lsh(0)[0]*j];
            }
            
            /// same as above but using iterator
            T operator()(const const_iter& i) const 
            {
               return (*this)(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// same as above but linear index
            T operator()(const int p, const int ij) const
            {
               assert(p < npatches);
               return data[ij];
            }
            
            /// modify local surface data on patch "p", index i,j
            T& operator()(const int p, const int i, const int j)
            {
               assert(p < npatches); 
               return data[i + lsh(0)[0]*j];
            }
            
            /// same as above but using iterator
            T& operator()(const const_iter& i)
            {
               return (*this)(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// same as above but linear index
            T& operator()(const int p, const int ij)
            {
               assert(p < npatches);
               return data[ij];
            }
            
            
            /// return pointer to data
            void* data_pointer()
            {
               return (void*) &data.front(); 
            }
            
            /// access local surface radius on patch "p", index i,j
            CCTK_REAL radius(const int p=0, const int i=0, const int j=0) const
            {
               assert(p < npatches); 
               if (spheredata<T>::has_constant_radius())
                  return spheredata<T>::_radius;
               else
               {
                  assert(_radii != NULL);
                  return _radii[i+lbnd(p)[0] + gsh(0)[0]*(j+lbnd(p)[1])]; // the radius is always stored completely on each processor!
               }
            }
            
            /// same as above: access local surface radius using iterator
            CCTK_REAL radius(const const_iter& i) const 
            {
               return radius(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// modify local surface radius on patch "p", index i,j
            CCTK_REAL& radius(const int p=0, const int i=0, const int j=0)
            {
               assert(p < npatches); 
               if (spheredata<T>::has_constant_radius())
                  return spheredata<T>::_radius;
               else
               {
                  assert(_radii != NULL);
                  return _radii[i+lbnd(p)[0] + gsh(0)[0]*(j+lbnd(p)[1])]; // the radius is always stored completely on each processor!
               }
            }
            
            /// same as above: modify local surface radius using iterator
            CCTK_REAL& radius(const const_iter& i)
            {
               return radius(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// return pointer to surface radius data
            void* radius_pointer() const 
            {
               // if this slice carries the radius-function
               // then return the data!
               if (spheredata<T>::varname() == "ss_radius")
                  return (void*) &data.front(); 
               return (void*)_radii;
            }
            
            /// access delta-spacing
            vect<CCTK_REAL, 2> delta() const { return vect<CCTK_REAL,2>(_dtheta, _dphi); }
            
            /// access local angular coordinates on patch p
            vect<CCTK_REAL, 2> coord(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
               return vect<CCTK_REAL,2>(_origin_theta+(i+_lbnd[0])*_dtheta, _origin_phi+(j+_lbnd[1])*_dphi); 
            }
            
            /// same as above but using iterator
            vect<CCTK_REAL, 2> coord(const const_iter& i) const 
            {
               return coord(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// same as above but using linear index
            vect<CCTK_REAL, 2> coord(const int p, const int ij) const 
            {
               assert(p < npatches);
               const int i = ij % _lsh[0];
               const int j = ij / _lsh[0];
               return vect<CCTK_REAL,2>(_origin_theta+(i+_lbnd[0])*_dtheta, _origin_phi+(j+_lbnd[1])*_dphi); 
            }
            
            /// access global angular coordinates (standard theta,phi spherical coordinate system which for this
            /// system is the same as the local coordinate system since we are already in theta,phi-coordinates)
            vect<CCTK_REAL, 2> coord_spherical(const int p, const int i, const int j) const { return coord(p, i, j); }
            
            /// same as above but using iterator
            vect<CCTK_REAL, 2> coord_spherical(const const_iter& i) const 
            {
               return coord_spherical(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns x-coordinate value of local point i,j on patch p
            CCTK_REAL cart_x(const int p, const int i, const int j) const 
            { 
               assert(p < npatches);
               return spheredata<T>::origin()[0] + radius(p, i, j) * sin(_origin_theta+(i+_lbnd[0])*_dtheta) * cos(_origin_phi+(j+_lbnd[1])*_dphi);
            }
            
            /// same as above but using iterator
            CCTK_REAL cart_x(const const_iter& i) const 
            {
               return cart_x(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns y-coordinate value of local point i,j on patch p
            CCTK_REAL cart_y(const int p, const int i, const int j) const 
            { 
               assert(p < npatches); 
               return spheredata<T>::origin()[1] + radius(p, i, j) * sin(_origin_theta+(i+_lbnd[0])*_dtheta) * sin(_origin_phi+(j+_lbnd[1])*_dphi);
            }
            
            /// same as above but using iterator
            CCTK_REAL cart_y(const const_iter& i) const 
            {
               return cart_y(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns z-coordinate value of local point i,j on patch p
            CCTK_REAL cart_z(const int p, const int i, const int j) const 
            { 
               assert(p < npatches); 
               return spheredata<T>::origin()[2] + radius(p, i, j) * cos(_origin_theta+(i+_lbnd[0])*_dtheta);
            }
            
            /// same as above but using iterator
            CCTK_REAL cart_z(const const_iter& i) const 
            {
               return cart_z(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// global surface size on patch "p"
            vect<int, 2> gsh(const int p) const { assert(p < npatches); return _gsh; }
            /// local size on patch "p"
            vect<int, 2> lsh(const int p) const { assert(p < npatches); return _lsh; }
            /// upper local bound on patch "p"
            vect<int, 2> ubnd(const int p) const { assert(p < npatches); return _ubnd; }
            /// upper local bound on patch "p"
            vect<int, 2> lbnd(const int p) const { assert(p < npatches); return _lbnd; }
            
            /// query whether given set of indices is in ghostzone
            bool ghostzone(const int p, const int i, const int j) const
            {
               if (i < spheredata<T>::nghosts() || 
                   i > lsh(p)[0]-spheredata<T>::nghosts()-1 ||
                   j < spheredata<T>::nghosts() || 
                   j > lsh(p)[1]-spheredata<T>::nghosts()-1)
                   return true;
               return false;
            }
            
            /// access any theta/phi coordinate via interpolation
            CCTK_REAL interpolate(const CCTK_REAL theta, const CCTK_REAL phi) const 
            { 
               // fourth-order interpolation
               /*const int order = 4
               const int n = order+1  // number of points
               int istart = 0
            
               // get point that is closest to our interpolation point
               for (int ii=0; ii < lsh(0)[0]; ++ii)
               {
                  if (coord(0, ii, 0)[0] <= theta and coord(0, ii+1, 0)[0] > theta)
                  {
                        if (fabs(coords[ii]-x) < fabs(coords[ii+1]-x))
                           istart = (int) (ii-round(n*1.0/2.0))
                        else
                           istart = (int) (ii+1-round(n*1.0/2.0))
                        break;
                  }
               }
               
               istart = (int) (ii-round(n*1.0/2.0))
               if (istart < 0) istart = 0
               if (istart+n >= lsh(0)[0]) istart = lsh(0)[0]-1 - n
               
               jstart = (int) (ii-round(n*1.0/2.0))
               if (jstart < 0) jstart = 0
               if (jstart+n >= lsh(0)[1]) jstart = lsh(0)[1]-1 - n
               
               val = 0
               
               for jj in range(istart, istart+n+1):
                  c_j = 1
                  for kk in range(istart, istart+n+1):
                        if (jj != kk):
                           c_j = c_j * (x - coords[kk]) / (coords[jj] - coords[kk])
                  val = val + c_j * data[jj]
               */
               return 0; 
            }
            
            /// interpolate from Cactus gridfunctions onto sphere
            void interpolate(const cGH* const cctkGH) 
            {
              // get the input variable's index
              const char* const varname = spheredata<T>::varname().c_str();
              vector<int> varindices (N_INPUT_ARRAYS, CCTK_VarIndex(varname));
              if (varindices[0] < 0) {
                CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                           "couldn't get index of slice variable '%s'", varname);
              }
              vector<CCTK_REAL *> values (N_OUTPUT_ARRAYS, &data.front());

              // set up the interpolation
              assert(interp_setups.size() > spheredata<T>::ID());
              interp_setup_t* &interp_setup = interp_setups[spheredata<T>::ID()];
              if (not interp_setup) {
                interp_setup = new interp_setup_t(lsh(0)[0] * lsh(0)[1]);

                // allocate storage for coordinates
                fasterp_glocs_t locations (interp_setup->npoints);

                // get Cartesian coordinate values of gridpoints on spherical surface
                for (int j=0; j < lsh(0)[1]; ++j) {
                  for (int k=0; k < lsh(0)[0]; ++k) {
                    const int l = k + lsh(0)[0]*j;
                    locations.coords[0][l] = cart_x(0, k, j);
                    locations.coords[1][l] = cart_y(0, k, j);
                    locations.coords[2][l] = cart_z(0, k, j);
                  }
                }

                interp_setup->fasterp_setup =
                  new fasterp_setup_t(cctkGH, locations, interpolator_order);
              }

              // do the interpolation
              assert(interp_setup->fasterp_setup);
              assert(interp_setup->npoints == lsh(0)[0] * lsh(0)[1]);
              interp_setup->fasterp_setup->interpolate (cctkGH, varindices, values);

              if (not spheredata<T>::has_constant_radius()) {
                delete interp_setup;
                interp_setup = NULL;
              }
            }
            
            /// surface integral over slice
            CCTK_REAL integrate() const 
            {
               CCTK_REAL result = 0;
               
               //if (data.size() > 0)
               {
                  // We do a Gauss-quadrature which is of order O(h^{2N})
                  // where N is the number of points.
                  if (spheredata<T>::has_constant_radius())
                  {
                     integrator myint(*this);
                     
                     for (const_iter i=begin(); !i.done(); ++i)
                        myint.sum(i, radius(0,0,0)*radius(0,0,0)*sin(i.idx().theta));
                     
                     result = myint.finalize();
                  }
                  else
                  {
                     // we need a ghostzone width of 2 for integration since we need to take
                     // derivatives!
                     if (spheredata<T>::nghosts() < 2) 
                        CCTK_WARN(0, "ghostzone-width too small!");
                  
                     integrator myint(*this);
                     
                     for (const_iter i=begin(); !i.done(); ++i)
                     {
                        // dont't integrate in ghostzones
                        if (i.ghostzone())
                           continue;
                     
                        myint.sum(i, det(i));
                     }
                     result = myint.finalize();
                  }
               }
               
               return result; 
            }

            
            /// take pointwise derivative on patch p and point i,j in theta direction
            CCTK_REAL dx(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
            
               if (i < 2 || i >= lsh(p)[0]-2) 
               {
                  CCTK_WARN(0, "ghostzone-width too small!");
                  return 0.0;
               }
               return ((1.0/(12.0*_dtheta)))*(-data[SINDEX2D(p, i+2, j)]
                                             + data[SINDEX2D(p, i-2, j)]
                                        + 8.0*(data[SINDEX2D(p, i+1, j)]
                                             - data[SINDEX2D(p, i-1, j)]));
            }
            
            /// same as above but using iterator
            CCTK_REAL dx(const const_iter& i) const 
            {
               return dx(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// take pointwise derivative on patch p and point i,j in phi direction
            CCTK_REAL dy(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
            
               if (j < 2 || j >= lsh(p)[1]-2) 
               {
                  CCTK_WARN(0, "ghostzone-width too small!");
                  return 0.0;
               }
               return ((1.0/(12.0*_dphi)))*(-data[SINDEX2D(p, i, j+2)]
                                           + data[SINDEX2D(p, i, j-2)]
                                      + 8.0*(data[SINDEX2D(p, i, j+1)]
                                           - data[SINDEX2D(p, i, j-1)]));
            }
            
            /// same as above but using iterator
            CCTK_REAL dy(const const_iter& i) const 
            {
               return dy(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dxdx(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
            
               if (i < 2 || i >= lsh(p)[0]-2) 
               {
                  CCTK_WARN(0, "ghostzone-width too small!");
                  return 0.0;
               }
               return ((1.0/(12.0*_dtheta*_dtheta))*(-(data[SINDEX2D(p, i+2, j)]
                                                     + data[SINDEX2D(p, i-2, j)])
                                               + 16.0*(data[SINDEX2D(p, i+1, j)]
                                                     + data[SINDEX2D(p, i-1, j)])
                                                     - 30.0*data[SINDEX2D(p, i, j)]));
            }
            
            /// same as above but using iterator
            CCTK_REAL dxdx(const const_iter& i) const 
            {
               return dxdx(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dxdy(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
            
               if (i < 2 || i >= lsh(p)[0]-2 ||
                   j < 2 || j >= lsh(p)[1]-2) 
               {
                  CCTK_WARN(0, "ghostzone-width too small!");
                  return 0.0;
               }
               return ((1.0/(144.0*_dtheta*_dphi)))*( (data[SINDEX2D(p, i+2, j+2)]
                                                     - data[SINDEX2D(p, i+2, j-2)]
                                                     - data[SINDEX2D(p, i-2, j+2)]
                                                     + data[SINDEX2D(p, i-2, j-2)])
                                               + 8.0*(-data[SINDEX2D(p, i+2, j+1)]
                                                      +data[SINDEX2D(p, i+2, j-1)]
                                                      -data[SINDEX2D(p, i+1, j+2)]
                                                      +data[SINDEX2D(p, i+1, j-2)]
                                                      -data[SINDEX2D(p, i-1, j-2)]
                                                      +data[SINDEX2D(p, i-1, j+2)]
                                                      -data[SINDEX2D(p, i-2, j-1)]
                                                      +data[SINDEX2D(p, i-2, j+1)])
                                               + 64.0*(data[SINDEX2D(p, i+1, j+1)]
                                                     - data[SINDEX2D(p, i+1, j-1)]
                                                     - data[SINDEX2D(p, i-1, j+1)]
                                                     + data[SINDEX2D(p, i-1, j-1)]));
            }
            
            /// same as above but using iterator
            CCTK_REAL dxdy(const const_iter& i) const 
            {
               return dxdy(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// take pointwise derivative on patch p and point i,j
            CCTK_REAL dydy(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
            
               if (j < 2 || j >= lsh(p)[1]-2) 
               {
                  CCTK_WARN(0, "ghostzone-width too small!");
                  return 0.0;
               }
               return ((1.0/(12.0*_dphi*_dphi))*(-(data[SINDEX2D(p, i, j+2)]
                                                 + data[SINDEX2D(p, i, j-2)])
                                           + 16.0*(data[SINDEX2D(p, i, j+1)]
                                                 + data[SINDEX2D(p, i, j-1)])
                                                 - 30.0*data[SINDEX2D(p, i, j)]));
            }
            
            /// same as above but using iterator
            CCTK_REAL dydy(const const_iter& i) const 
            {
               return dydy(i.idx().p, i.idx().i, i.idx().j);
            }
            
            
            /// returns the L2-norm over the sphere
            CCTK_REAL normL2() const
            {
               CCTK_REAL norm = 0.0;
               for (const_iter i=begin(); !i.done(); ++i)
                  if (!i.ghostzone())
                     norm += (*i) * (*i);
               // interprocessor synchronization
               CCTK_REAL global_norm;
               MPI_Allreduce (&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
               return sqrt(global_norm)/((spheredata<T>::gsh(0)[0]-2*spheredata<T>::nghosts())*(spheredata<T>::gsh(0)[1]-2*spheredata<T>::nghosts()));
            }

            /// returns the infinity-norm over the sphere
            CCTK_REAL normLinf() const
            {
               CCTK_REAL norm = *(begin());
               for (const_iter i=begin(); !i.done(); ++i)
                  if (fabs(*i) > norm)
                     norm = fabs(*i);
               // interprocessor synchronization
               CCTK_REAL global_norm;
               MPI_Allreduce (&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
               return global_norm;
            }
            
            
            /// project onto spin-weighted spherical harmonic with spin
            /// s, and usual l, m
            CCTK_COMPLEX contract_with_sYlm(const int s, const int l, const int m) const 
            {
               CCTK_REAL result_re = 0;
               CCTK_REAL result_im = 0;
               
               //if (data.size() > 0)
               {
                  // We do a Gauss-quadrature which is of order O(h^{2N})
                  // where N is the number of points.
                  if (spheredata<T>::has_constant_radius())
                  {
                     integrator intsYlm_re(*this);
                     integrator intsYlm_im(*this);
                     
                     for (const_iter i=begin(); !i.done(); ++i)
                     {
                        double sYlm_re, sYlm_im;
	                sYlm(s,l,m, i.idx().theta, i.idx().phi, &sYlm_re, &sYlm_im);
                        intsYlm_re.sum(i, sin(i.idx().theta), sYlm_re);
                        intsYlm_im.sum(i, sin(i.idx().theta), sYlm_im);
                     }
                     
                     result_re = intsYlm_re.finalize();
                     result_im = intsYlm_im.finalize();
                  }
                  else
                  {
                     // we need a ghostzone width of 2 for integration since we need to take
                     // derivatives!
                     if (spheredata<T>::nghosts() < 2) 
                        CCTK_WARN(0, "ghostzone-width too small!");
                  
                     integrator intsYlm_re(*this);
                     integrator intsYlm_im(*this);
                     
                     for (const_iter i=begin(); !i.done(); ++i)
                     {
                        // dont't integrate in ghostzones
                        if (i.ghostzone())
                           continue;
                        
                        CCTK_REAL r = radius(i.idx().p, i.idx().i, i.idx().j);
                        double sYlm_re, sYlm_im;
	                sYlm(s,l,m, i.idx().theta, i.idx().phi, &sYlm_re, &sYlm_im);
                        intsYlm_re.sum(i, det(i) / (r*r), sYlm_re);
                        intsYlm_im.sum(i, det(i) / (r*r), sYlm_im);
                     }
                     
                     result_re = intsYlm_re.finalize();
                     result_im = intsYlm_im.finalize();
                  }
               }
               
               return CCTK_Cmplx(result_re, -result_im); 
            }
            
            
            /// calculate surface determinant for given point
            CCTK_REAL det(const int p, const int i, const int j) const
            {
               assert(p < npatches);
               
               const CCTK_REAL theta = _origin_theta+(i+_lbnd[0])*_dtheta;
               
               if (spheredata<T>::has_constant_radius())
                  return radius(0,0,0)*radius(0,0,0)*sin(theta);
                  
               assert(i >= 2 && i < lsh(p)[0]-2);
               assert(j >= 2 && j < lsh(p)[1]-2);
               
               const CCTK_REAL r = radius(p, i, j);
                        
               const CCTK_REAL drdtheta = ((1.0/(12.0*_dtheta)))*(-radius(p, i+2, j)
                                                           + radius(p, i-2, j)
                                                      + 8.0*(radius(p, i+1, j)
                                                           - radius(p, i-1, j)));
               
               
               const CCTK_REAL drdphi = ((1.0/(12.0*_dphi)))*(-radius(p, i, j+2)
                                                       + radius(p, i, j-2)
                                                  + 8.0*(radius(p, i, j+1)
                                                       - radius(p, i, j-1)));
               
               
               return sqrt(r*r*(drdphi*drdphi + sin(theta)*sin(theta)*(r*r + drdtheta*drdtheta) ));
            }
            
            /// same as above but with iterator
            CCTK_REAL det(const_iter& i) const
            {
               return det(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// write the SphereData including all attributes
            /// to the output stream defined by io_base (or its inheritants)
            void output (io_base& io) const ;
            
            void input (io_base& io) { }
            
            iter begin() { return iter(this); }
            const_iter begin() const { return const_iter(this); }
            
            
            /// iterator class to traverse the grid of the slice
            class const_iter
            {
               public :
                        const_iter(spheredata_1patch const * const spheredata_) 
                           : _spheredata(spheredata_), 
                              _idx(0,0,0,0,spheredata_->_origin_theta+ (spheredata_->lbnd(0)[0]*spheredata_->_dtheta), 
                                           spheredata_->_origin_phi + (spheredata_->lbnd(0)[1]*spheredata_->_dphi)) { }
                        virtual ~const_iter() { }
                        
                        /// dereferencing operator as data-point accessor
                        T operator*() const 
                        {
                           assert(_spheredata->data.size() > 0); 
                           return (*_spheredata).data[_idx.ij];  
                        }
                        
                        spheredata_1patch<T> const * spheredata() const
                        {
                           return _spheredata; 
                        }
                        
                        void operator++() { (*this)++; }
                        
                        /// increment iterator
                        void operator++(int)
                        {
                           if (_idx.i >= _spheredata->lsh(_idx.p)[0]-1) 
                           {
                              _idx.i=0;
                              _idx.theta = _spheredata->_origin_theta + (_spheredata->lbnd(_idx.p)[0]*_spheredata->_dtheta);
                              if (_idx.j >= _spheredata->lsh(_idx.p)[1]-1)
                              {
                                 _idx.p++;
                                 if (_idx.p == _spheredata->npatches)
                                    return;
                                 _idx.theta = _spheredata->_origin_theta + (_spheredata->lbnd(_idx.p)[0]*_spheredata->_dtheta);
                                 _idx.phi = _spheredata->_origin_phi + (_spheredata->lbnd(_idx.p)[1]*_spheredata->_dphi);
                                 _idx.ij = 0;
                                 _idx.i = 0;
                                 _idx.j = 0;
                              }
                              else
                              {
                                 _idx.j++;
                                 _idx.ij++;
                                 _idx.phi+=_spheredata->_dphi;
                              }
                           }
                           else
                           {
                              _idx.i++;
                              _idx.ij++;
                              _idx.theta+=_spheredata->_dtheta;
                           }
                        }
                        
                        /// query whether iterator is done
                        bool done() const
                        {
                           if (_idx.p == _spheredata->npatches)
                              return true;
                           return false;
                        }
                        
                        /// query whether iterator is in ghostzone
                        bool ghostzone() const
                        {
                           /*if (_idx.i < _spheredata->nghosts() || 
                               _idx.i > _spheredata->lsh(_idx.p)[0]-_spheredata->nghosts()-1 ||
                               _idx.j < _spheredata->nghosts() || 
                               _idx.j > _spheredata->lsh(_idx.p)[1]-_spheredata->nghosts()-1)
                               return true;
                           return false;*/
                           return _spheredata->ghostzone(_idx.p, _idx.i, _idx.j);
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
                           
                           idx_t(const int p_, 
                                 const int i_, 
                                 const int j_, 
                                 const int ij_, 
                                 const CCTK_REAL theta_, 
                                 const CCTK_REAL phi_) : p(p_), i(i_), j(j_), ij(ij_), theta(theta_), phi(phi_) { }
                        };
                        
                        /// access current index-struct
                        idx_t idx() const { return _idx; }
                        
               protected :
                        idx_t _idx;
               private :
                        spheredata_1patch<T> const * const _spheredata;
            };
            
            
            class iter : public const_iter 
            {
               public :
                        iter(spheredata_1patch* const spheredata_) : _spheredata(spheredata_), const_iter(spheredata_) { }
                        virtual ~iter() { }
                        
                        /// dereferencing operator as data-point accessor
                        T& operator*() 
                        {
                           assert(_spheredata->data.size() > 0);  
                           return (*_spheredata).data[this->_idx.ij];  
                        }
                        
                        spheredata_1patch<T>* spheredata()
                        {
                           return _spheredata; 
                        }
               private :
                        spheredata_1patch<T>* const _spheredata;
            };
            
            
            /// integrator class
            class integrator
            {
               public :
                        integrator(const spheredata_1patch<T>& sph_) : sph(sph_) { init(); }
                        virtual ~integrator() { }
                        
                        /// initialize
                        void init() 
                        {
                            theta_sum = 0.0;
                            _result = 0.0;
                            finalized = false;
                        }
                        /// finalize integration and return result
                        CCTK_REAL finalize(commstack* const cs = NULL) 
                        {
                           finalized = true;
                           _result *= sph.delta()[0]*8.0/(sph.gsh(0)[1]-2*sph.nghosts());
                           
                           // synchronize results!
                           if (!cs)
                           {
                              if (sph.distrib_method() != constant)
                              {
                                 double dummy; 
                                 MPI_Allreduce (&_result, &dummy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                                 _result = dummy;
                              }
                           }
                           else
                           {
                              cs->push(_result);
                              return -1;
                           }
                           
                           
                           return _result;
                        }
                        
                        /// return the result of the integration (integration must already be finalized)
                        CCTK_REAL result() const
                        {
                           assert(finalized == true);
                           return _result;
                        }
                        
                        /// sum over function on the sphere (and multiply by given function and determinant)
                        void sum(const int p, const int i, const int j, const CCTK_REAL f = 1.0, const CCTK_REAL det = 1.0) 
                        {
                           // dont't integrate in ghostzones
                           if (sph.ghostzone(p,i,j))
                              return;
                           
                           if (i == sph.nghosts())
                              theta_sum = 0;
                           
                           theta_sum += f * det * sph.int_weights[i-sph.nghosts()];
                                     //* weight(it.idx().i+sph.lbnd(0)[0]-sph.nghosts(), sph.gsh(0)[0]-2*sph.nghosts());
                           
                           if (i == sph.lsh(p)[0]-sph.nghosts()-1)
                              _result += theta_sum;
                        }
                        
                        /// same as above except that there is also an additive function now
                        void sum(const int p, const int i, const int j, const CCTK_REAL f, const CCTK_REAL g, const CCTK_REAL det) 
                        {
                           // dont't integrate in ghostzones
                           if (sph.ghostzone(p,i,j))
                              return;
                           
                           if (i == sph.nghosts())
                              theta_sum = 0;
                           
                           theta_sum += (f + g) * det * sph.int_weights[i-sph.nghosts()];
                                     //* weight(it.idx().i+sph.lbnd(0)[0]-sph.nghosts(), sph.gsh(0)[0]-2*sph.nghosts());
                           
                           if (i == sph.lsh(p)[0]-sph.nghosts()-1)
                              _result += theta_sum;
                           
                           // Why not simply doing it like this???
                           //_result += (f + g) * det * sph.int_weights[i-sph.nghosts()];
                        }
                        
                        /// sum over function on the sphere (and multiply by given function and determinant)
                        void sum(const_iter& it, const CCTK_REAL f = 1.0, const CCTK_REAL det = 1.0) 
                        {
                           return sum(it.idx().p, it.idx().i, it.idx().j, (*it)*f, det);
                        }
                        
                        
                        void sum(const_iter& it, const CCTK_REAL f, const CCTK_REAL g, const CCTK_REAL det)
                        {
                           return sum(it.idx().p, it.idx().i, it.idx().j, (*it)*f, g, det);
                        }
                        
                        
               private :
                        const spheredata_1patch<T>& sph;
                        CCTK_REAL theta_sum;
                        CCTK_REAL _result;
                        bool finalized;
            };
            
   private :
            
            /// decompose local piece of data for this MPI process
            void decompose()
            {
               // we need to make sure that slices belonging to the same
               // slice-no. in the parfile are split in exactly the same way!
               // Otherwise simultanoues gridpoint loops with many slices of the same slice-no. may get spoiled!
               
               int nprocs = 1;
               MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
               
               // currently, we will split accross all processors!
               assert(nprocs == spheredata<T>::_processors.size());
               
               vect<int,2> dims(0, 0);
               MPI_Dims_create(nprocs, 2, &dims[0]);
               int used_procs = dims[0]*dims[1];
               
               //cout << "-------- dims = " << dims << endl; 
               
               int this_proc = 0;
               MPI_Comm_rank(MPI_COMM_WORLD, &this_proc);
               
               vect<int,2> proc_vect = multi_index<2>(this_proc, dims);
               
               //cout << "-------- proc_vect = " << proc_vect << endl; 
               
               if (this_proc < used_procs)
               {
                  _lsh = vect<int,2>(0, 0);
                  _lbnd = vect<int,2>(0, 0);
                  _ubnd = vect<int,2>(0, 0);
                  
                  _lsh = (_gsh - vect<int,2>(2*spheredata<T>::nghosts(), 2*spheredata<T>::nghosts())) / dims + 2*spheredata<T>::nghosts();
                  
                  vect<int,2> mod = (_gsh - vect<int,2>(2*spheredata<T>::nghosts(), 2*spheredata<T>::nghosts())) % dims;
                  
                  for (int i=0; i < 2; ++i)
                     if (proc_vect[i] < mod[i])
                        _lsh[i]++;
                  
                  assert(_lsh[0] != 2*spheredata<T>::nghosts() && _lsh[1] != 2*spheredata<T>::nghosts());
                  
                  for (int i=0; i < 2; i++)
                  {
                     if (proc_vect[i] == 0)
                        _lbnd[i] = 0;
                     else
                     {
                        _lbnd[i] = _lsh[i] * proc_vect[i] - proc_vect[i]*2*spheredata<T>::nghosts();
                        if (proc_vect[i] >= mod[i])
                           _lbnd[i] += mod[i];
                     }
                     if (proc_vect[i] == dims[i]-1)
                        _ubnd[i] = _gsh[i]-1;
                     else
                        _ubnd[i] = _lbnd[i] + _lsh[i]-1;// + spheredata<T>::nghosts();
                  }
      
                  /*for (int i=0; i < 2; i++)
                  {
                     if (proc_vect[i] == dims[i]-1)
                     {
                        _lsh[i] += (_gsh[i]- 2*spheredata<T>::nghosts()) % dims[i];
                        _ubnd[i] += (_gsh[i]- 2*spheredata<T>::nghosts()) % dims[i];
                     }
                  }*/
                  
                  //cout << "----------- " << this_proc << ": " << _lsh << ",  " << _lbnd << ",   " << _ubnd << endl;
                  
                  return;
               }
               
               // this processor carries no data....
               _lsh = vect<int,2>(0, 0);
               _lbnd = vect<int,2>(0, 0);
               _ubnd = vect<int,2>(0, 0);
            }
            
            
            /// the data of the entire sphere
            vector<T> data;
            
            /// the pointer to the radii of this slice number at each point if requested
            CCTK_REAL* _radii;
            
            /// surface descriptors
            vect<int, 2> _gsh;
            vect<int, 2> _lsh;
            vect<int, 2> _lbnd;
            vect<int, 2> _ubnd;
            
            /// the spacing on the sphere
            CCTK_REAL _dphi, _dtheta;
            
            /// the angular origin
            CCTK_REAL _origin_theta, _origin_phi;
            
            /// precalculated weights for Gauss-integration along theta direction (1d array).
            /// Without this, integration becomes a real bottleneck!
            vector<CCTK_REAL> int_weights;
};



template <typename T>
inline io_base& operator<< (io_base& io, const spheredata_1patch<T>& sd)
{
   sd.output(io);
   
   return io;
} 



}


#endif

