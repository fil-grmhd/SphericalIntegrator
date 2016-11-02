
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

#ifndef _SPHEREDATA6_
#define _SPHEREDATA6_

#include <vector>

/// use Carpet-vectors for fast and easy vector-handling
#include "vect.hh"

#include "spheredata.hh"
#include "commstack.hh"
#include <math.h>


namespace SPS {

using namespace std;


#define PLUS_X  0
#define MINUS_X 1
#define PLUS_Y  2
#define MINUS_Y 3
#define PLUS_Z  4
#define MINUS_Z 5


/**
   A 6-patch system covering the sphere.
   It offers a uniform spatial sampling of the sphere.
   In addition it will work hand in hand with Llama
   (Thornburg04 coordinates).
*/
template <typename T>
class spheredata_6patch : public spheredata<T>
{
   public :
            static int const npatches = 6;
            class iter;  // forward declaration
            class const_iter;  // forward declaration
            
            spheredata_6patch() : _desc(0) { }
            spheredata_6patch(const spheredata_6patch& sd) 
            { 
               *this = sd;
               if (spheredata<T>::varname() == "ss_radius")
                  _radii = &data; 
            }
            spheredata_6patch(const string& varname_,
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
               : _desc(0),
                 _radii((vector<vector<CCTK_REAL> >*) radii_),
                 spheredata<T>(varname_, id_, ntheta_, nphi_, nghosts_, radius_, radii_,
                               origin_, has_constant_radius_, symmetry_,
                               distrib_method_, processors_, can_use_Llama_)
            {
               _desc.resize(npatches);
               data.resize(npatches);
               //_radii.resize(npatches);
               
               if (spheredata<T>::_distrib_method == constant) // all processors will carry all data
               {
                  for (int i=0; i < npatches; ++i)
                  {
                     _desc[i] = patch_desc(spheredata<T>::npoints(), 
                                           spheredata<T>::npoints(),
                                           vect<int,2>(0),
                                           vect<int,2>(spheredata<T>::ntheta(), spheredata<T>::nphi()),
                                           spheredata<T>::nghosts());
                  
                     data[i].resize(_desc[i].lsh()[0]*_desc[i].lsh()[1]);
                  }
               }
               if (spheredata<T>::_distrib_method == single) // one processor will carry all data
               {
                  int this_proc = 0;
                  assert(spheredata<T>::_processors.size() == 1);
                  MPI_Comm_rank(MPI_COMM_WORLD, &this_proc);
                  if (this_proc == spheredata<T>::_processors[0])
                  {
                     for (int i=0; i < npatches; ++i)
                     {
                        _desc[i] = patch_desc(spheredata<T>::npoints(), 
                                              spheredata<T>::npoints(),
                                              vect<int,2>(0),
                                              vect<int,2>(spheredata<T>::ntheta(), spheredata<T>::nphi()),
                                              spheredata<T>::nghosts());
                     
                        data[i].resize(_desc[i].lsh()[0]*_desc[i].lsh()[1]);
                     }
                  }
                  else
                  {
                     for (int i=0; i < npatches; ++i)
                     {
                        _desc[i] = patch_desc(spheredata<T>::npoints(), 
                                              0,
                                              0,
                                              0,
                                              0);
                        
                        data[i].resize(0);
                     }
                  }
               }
               if (spheredata<T>::_distrib_method == split) // a group of processors will carry the data
               {
                  for (int i=0; i < npatches; ++i)
                  {
                     _desc[i] = patch_desc(spheredata<T>::npoints(), 
                                           spheredata<T>::npoints(),
                                           vect<int,2>(0),
                                           vect<int,2>(spheredata<T>::ntheta(), spheredata<T>::nphi()),
                                           spheredata<T>::nghosts());
                  
                  }
               
                  // do a domain decomposition (modifies patch_desc)
                  decompose();
                  
                  for (int i=0; i < npatches; ++i)
                     data[i].resize(_desc[i].lsh()[0]*_desc[i].lsh()[1]);
               }
               
               // set to poison values initially
               for (int j=0; j < npatches; ++j)
                  for (int i=0; i < data[j].size(); ++i)
                     data[j][i] = POISON_VAL;
               
               // if the radius-pointer has not been set then we assume that this slice
               // cariies the radius as data!
               if (_radii == NULL)
               {
                  _radii = &data;
                  
                  // set radius according to initial constant value
                  for (iter i=begin(); !i.done(); ++i)
                     *i = radius_;
               }
            }
            
            virtual ~spheredata_6patch() { }
            
            /// determines whether this slice takes advantage of Llama
            /// (only the 6patch slices can potentially take advantage)
            bool uses_Llama() const { return _uses_Llama; }
            
            /// given two surface indices this will return the linear index
            int SINDEX2D(const int p, const int i, const int j) const
            {
               return i + lsh(p)[0]*j;
            }
            
            /// access local surface data on patch "p", index i,j
            T operator()(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
               return data[p][i + lsh(p)[0]*j];
            }
            
            /// same as above but using iterator
            T operator()(const const_iter& i) const 
            {
               return (*this)(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// modify local surface data on patch "p", index i,j
            T& operator()(const int p, const int i, const int j)
            {
               assert(p < npatches); 
               return data[p][i + lsh(p)[0]*j];
            }
            
            /// same as above but using iterator
            T& operator()(const const_iter& i)
            {
               return (*this)(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// return pointer to data
            void* data_pointer(const int m)
            {
               return (void*) &data[m].front(); 
            }
            
            /// access local surface radius on patch "p", index i,j
            CCTK_REAL radius(const int p, const int i, const int j) const
            {
               assert(p < npatches); 
               if (spheredata<T>::has_constant_radius())
                  return spheredata<T>::_radius;
               else
               {
                  assert(_radii != NULL);
                  return (*_radii)[p][i+lbnd(p)[0] + gsh(0)[0]*(j+lbnd(p)[1])]; // the radius is always stored completely on each processor!
               }
            }
            
            /// same as above: access local surface radius using iterator
            CCTK_REAL radius(const const_iter& i) const 
            {
               return radius(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// modify local surface radius on patch "p", index i,j
            CCTK_REAL& radius(const int p, const int i, const int j)
            {
               assert(p < npatches); 
               if (spheredata<T>::has_constant_radius())
                  return spheredata<T>::_radius;
               else
               {
                  assert(_radii != NULL);
                  return (*_radii)[p][i+lbnd(p)[0] + gsh(0)[0]*(j+lbnd(p)[1])]; // the radius is always stored completely on each processor!
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
               if (spheredata<T>::_varname == "ss_radius")
                  return (void*) &data; 
               return (void*)_radii; 
            }
            
            /// access delta-spacing
            vect<CCTK_REAL, 2> delta() const { return vect<CCTK_REAL,2>(_desc[0].dtheta(), _desc[0].dphi()); }
            
            /// access local angular coordinates on patch p (local 6-patch ("inflated cube") coordinates)
            vect<CCTK_REAL, 2> coord(const int p, const int i, const int j) const 
            {
               assert(p < npatches);
               return vect<CCTK_REAL,2>(_desc[p].origin_theta()+(i+lbnd(p)[0])*_desc[p].dtheta(), _desc[p].origin_phi()+(j+lbnd(p)[1])*_desc[p].dphi()); 
            }
            
            /// same as above but using iterator
            vect<CCTK_REAL, 2> coord(const const_iter& i) const 
            {
               return coord(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// access global angular coordinates (standard theta,phi spherical coordinate system)
            vect<CCTK_REAL, 2> coord_spherical(const int p, const int i, const int j) const;
            
            /// same as above but using iterator
            vect<CCTK_REAL, 2> coord_spherical(const const_iter& i) const 
            {
               return coord_spherical(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns x-coordinate value of local point i,j on patch p
            CCTK_REAL cart_x(const int p, const int i, const int j) const;
            
            /// same as above but using iterator
            CCTK_REAL cart_x(const const_iter& i) const 
            {
               return cart_x(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns y-coordinate value of local point i,j on patch p
            CCTK_REAL cart_y(const int p, const int i, const int j) const;
            
            /// same as above but using iterator
            CCTK_REAL cart_y(const const_iter& i) const 
            {
               return cart_y(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// returns z-coordinate value of local point i,j on patch p
            CCTK_REAL cart_z(const int p, const int i, const int j) const;
            
            /// same as above but using iterator
            CCTK_REAL cart_z(const const_iter& i) const 
            {
               return cart_z(i.idx().p, i.idx().i, i.idx().j);
            }
            
            /// global surface size on patch "p"
            vect<int, 2> gsh(const int p) const { assert(p < npatches); return _desc[p].gsh(); }
            /// local size on patch "p"
            vect<int, 2> lsh(const int p) const { assert(p < npatches); return _desc[p].lsh(); }
            /// upper local bound on patch "p"
            vect<int, 2> ubnd(const int p) const { assert(p < npatches); return _desc[p].ubnd(); }
            /// upper local bound on patch "p"
            vect<int, 2> lbnd(const int p) const { assert(p < npatches); return _desc[p].lbnd(); }
            
            /// access any theta/phi coordinate via interpolation
            CCTK_REAL interpolate(const CCTK_REAL theta, const CCTK_REAL phi) const { return 0; }
            
            /// interpolate from Cactus gridfunctions onto sphere
            void interpolate(const cGH* const cctkGH) 
            {
               const void* interp_coords[N_DIMS];
               CCTK_INT    input_array_indices[N_INPUT_ARRAYS];
               void*       output_arrays[N_OUTPUT_ARRAYS];
               const CCTK_INT input_array_type_codes[N_INPUT_ARRAYS]
                        = { CCTK_VARIABLE_REAL };
               const CCTK_INT output_array_type_codes[N_OUTPUT_ARRAYS]
                        = { CCTK_VARIABLE_REAL };
               
               // set interpolation handles
               int operator_handle = CCTK_InterpHandle("Lagrange polynomial interpolation");
               if (operator_handle < 0)
                  CCTK_WARN(0, "can't get interpolation handle!");
               int param_table_handle = Util_TableCreateFromString("order=4");
               if (param_table_handle < 0)
                  CCTK_WARN(0, "can't create parameter table!");
               int coordsys_handle = CCTK_CoordSystemHandle("cart3d");
               if (coordsys_handle < 0)
                  CCTK_WARN(0, "can't create coordsys handle! Forgot to activate a Coordinate-Thorn?");
               
               input_array_indices[0] = CCTK_VarIndex(spheredata<T>::varname().c_str());
               if (input_array_indices[0] < 0)
                  CCTK_WARN(0, "error getting VarIndex of variable that shall be sliced");
               
               for (int m=0; m < npatches; ++m)
               {
                  // project variable onto sphere
                  // this means we have to 3d-interpolate every non-coinciding gridpoint of
                  // the spherical grid from neighboring gridpoints of the Cartesian grid.
                  // for this, we first calculate the Cartesian coordinate of the jk-th spherical
                  // surface gridpoint
                  int n_interp_points = lsh(m)[0]*lsh(m)[1];
                  
                  // Arrays of Cartesian coordinates of the surface points onto which we want to interpolate
                  vector<CCTK_REAL> interp_x(n_interp_points, POISON_VAL);
                  vector<CCTK_REAL> interp_y(n_interp_points, POISON_VAL);
                  vector<CCTK_REAL> interp_z(n_interp_points, POISON_VAL);
                  
#ifdef DEBUG
                  CCTK_INFO("6patch: Getting Cartesian coords.");
#endif
                  // get Cartesian coordinate values of gridpoints on spherical surface
                  for (int j=0; j < lsh(0)[1]; ++j)
                     for (int k=0; k < lsh(0)[0]; ++k)
                     {
                        const int l = k + lsh(m)[0]*j;
                        interp_x[l] = cart_x(m, k, j);
                        interp_y[l] = cart_y(m, k, j);
                        interp_z[l] = cart_z(m, k, j);
                     }
                  
                  interp_coords[0] = (const void*) &interp_x.front();
                  interp_coords[1] = (const void*) &interp_y.front();
                  interp_coords[2] = (const void*) &interp_z.front();
                  
                  output_arrays[0] = (void*) &data[m].front();
            
#ifdef DEBUG
                  CCTK_INFO("6patch: Calling Interpolator");
#endif
                  // Do the actual interpolation.
                  // Only those processes interpolate that contain data
                  // All other processes simply send their
                  // grid-data.
                  if (data[m].size() == 0)
                     n_interp_points = 0;      // all other processes shall not interpolate
                  
                  if (CCTK_InterpGridArrays(cctkGH,
                                            N_DIMS,
                                            operator_handle, param_table_handle,
                                            coordsys_handle,
                                            n_interp_points,
                                               CCTK_VARIABLE_REAL,
                                               interp_coords,
                                            N_INPUT_ARRAYS,
                                               input_array_indices,
                                            N_OUTPUT_ARRAYS,
                                               output_array_type_codes,
                                               output_arrays) < 0)
                     CCTK_WARN(1, "error return from interpolator!");
                     
#ifdef DEBUG
                  CCTK_INFO("6patch: Done interpolating");
#endif
               }
            }
            
            /// surface integral over slice
            CCTK_REAL integrate() const 
            {
               CCTK_REAL result = 0;
               
               
               if (spheredata<T>::has_constant_radius())
               {
                  /*for (int j=spheredata<T>::nghosts(); j < lsh(m)[1]-spheredata<T>::nghosts(); ++j)
                  {
                     CCTK_REAL weight_phi = 1.0;
                     
                     if (j == spheredata<T>::nghosts() || j == lsh(m)[1]-spheredata<T>::nghosts()-1)
                        weight_phi = 0.5;
                     
                     CCTK_REAL theta_sum = 0;
                     for (int i=spheredata<T>::nghosts(); i < lsh(m)[0]-spheredata<T>::nghosts(); ++i)
                     {
                        CCTK_REAL weight_theta = 1.0;
                        
                        if (i == spheredata<T>::nghosts() || i == lsh(m)[0]-spheredata<T>::nghosts()-1)
                           weight_theta = 0.5;
                        
                        theta_sum += weight_theta * weight_phi * (*this)(m, i, j)  * det(m, i, j);
                     }
                     result_singlepatch += theta_sum;
                  }
                  result_singlepatch *= _desc[m].dtheta()*_desc[m].dphi()*radius(0,0,0)*radius(0,0,0);*/
                  
                  integrator myint(*this);
               
                  for (const_iter i=begin(); !i.done(); ++i)
                     myint.sum(i, det(i));
                  
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
               return ((1.0/(12.0*_desc[p].dtheta())))*(-data[p][SINDEX2D(p, i+2, j)]
                                             + data[p][SINDEX2D(p, i-2, j)]
                                        + 8.0*(data[p][SINDEX2D(p, i+1, j)]
                                             - data[p][SINDEX2D(p, i-1, j)]));
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
               return ((1.0/(12.0*_desc[p].dphi())))*(-data[p][SINDEX2D(p, i, j+2)]
                                           + data[p][SINDEX2D(p, i, j-2)]
                                      + 8.0*(data[p][SINDEX2D(p, i, j+1)]
                                           - data[p][SINDEX2D(p, i, j-1)]));
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
               return ((1.0/(12.0*_desc[p].dtheta()))*(-(data[p][SINDEX2D(p, i+2, j)]
                                                     + data[p][SINDEX2D(p, i-2, j)])
                                               + 16.0*(data[p][SINDEX2D(p, i+1, j)]
                                                     + data[p][SINDEX2D(p, i-1, j)])
                                                - 30.0*data[p][SINDEX2D(p, i, j)]));
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
               return ((1.0/(144.0*_desc[p].dtheta()*_desc[p].dphi())))*( (data[p][SINDEX2D(p, i+2, j+2)]
                                                     - data[p][SINDEX2D(p, i+2, j-2)]
                                                     - data[p][SINDEX2D(p, i-2, j+2)]
                                                     + data[p][SINDEX2D(p, i-2, j-2)])
                                               + 8.0*(-data[p][SINDEX2D(p, i+2, j+1)]
                                                      +data[p][SINDEX2D(p, i+2, j-1)]
                                                      -data[p][SINDEX2D(p, i+1, j+2)]
                                                      +data[p][SINDEX2D(p, i+1, j-2)]
                                                      -data[p][SINDEX2D(p, i-1, j-2)]
                                                      +data[p][SINDEX2D(p, i-1, j+2)]
                                                      -data[p][SINDEX2D(p, i-2, j-1)]
                                                      +data[p][SINDEX2D(p, i-2, j+1)])
                                               + 64.0*(data[p][SINDEX2D(p, i+1, j+1)]
                                                     - data[p][SINDEX2D(p, i+1, j-1)]
                                                     - data[p][SINDEX2D(p, i-1, j+1)]
                                                     + data[p][SINDEX2D(p, i-1, j-1)]));
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
               return ((1.0/(12.0*_desc[p].dphi()*_desc[p].dphi()))*(-(data[p][SINDEX2D(p, i, j+2)]
                                                 + data[p][SINDEX2D(p, i, j-2)])
                                           + 16.0*(data[p][SINDEX2D(p, i, j+1)]
                                                 + data[p][SINDEX2D(p, i, j-1)])
                                            - 30.0*data[p][SINDEX2D(p, i, j)]));
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
#warning "Don't know which patch index to use below -- using -1 instead"
               assert (0);
               return sqrt(global_norm)/(6*(spheredata<T>::gsh(-1)[0]-2*spheredata<T>::nghosts())*(spheredata<T>::gsh(-1)[1]-2*spheredata<T>::nghosts()));
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
               
               if (spheredata<T>::has_constant_radius())
               {
                  integrator intsYlm_re(*this);
                  integrator intsYlm_im(*this);
                  
                  for (const_iter i=begin(); !i.done(); ++i)
                  {
                     CCTK_REAL r = radius(i.idx().p, i.idx().i, i.idx().j);
                     double sYlm_re, sYlm_im;
	             sYlm(s,l,m, coord_spherical(i.idx().p, i.idx().i, i.idx().j)[0],
                                 coord_spherical(i.idx().p, i.idx().i, i.idx().j)[1],
			         &sYlm_re, &sYlm_im);
                     
                     intsYlm_re.sum(i, det(i) / (r*r), 
                                       sYlm_re);
                     intsYlm_im.sum(i, det(i) / (r*r), 
                                       sYlm_im);
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
	             sYlm(s,l,m, coord_spherical(i.idx().p, i.idx().i, i.idx().j)[0],
                                 coord_spherical(i.idx().p, i.idx().i, i.idx().j)[1],
			         &sYlm_re, &sYlm_im);
                     
                     intsYlm_re.sum(i, det(i) / (r*r), 
                                       sYlm_re);
                     intsYlm_im.sum(i, det(i) / (r*r), 
                                       sYlm_im);
                  }
                  
                  result_re = intsYlm_re.finalize();
                  result_im = intsYlm_im.finalize();
               }
               
               return CCTK_Cmplx(result_re, -result_im); 
            }
            
            
            /// calculates surface determinant for given point
            CCTK_REAL det(const_iter& i) const
            {
               assert(i.idx().p < npatches);
               
               CCTK_REAL rho = i.idx().theta;
               CCTK_REAL sigma = i.idx().phi;
               
               if (spheredata<T>::has_constant_radius())
                  return radius(0,0,0)*radius(0,0,0)*sqrt((pow(1.0/cos(rho),6) * pow(1.0/cos(sigma),6) - pow(1.0/cos(rho),4) * pow(1.0/cos(sigma),4) * pow(tan(rho),2) * pow(tan(sigma),2)))
                                                         / pow(pow(1.0/cos(sigma),2) + pow(tan(rho),2),2);
               
               assert(i.idx().i >= 2 && i.idx().i < lsh(i.idx().p)[0]-2);
               assert(i.idx().j >= 2 && i.idx().j < lsh(i.idx().p)[1]-2);
               
               CCTK_REAL r = radius(i.idx().p, i.idx().i, i.idx().j);
                        
               CCTK_REAL drdrho = ((1.0/(12.0*_desc[i.idx().p].dtheta())))*(-radius(i.idx().p, i.idx().i+2, i.idx().j)
                                                                   + radius(i.idx().p, i.idx().i-2, i.idx().j)
                                                              + 8.0*(radius(i.idx().p, i.idx().i+1, i.idx().j)
                                                                   - radius(i.idx().p, i.idx().i-1, i.idx().j)));
               
               
               CCTK_REAL drdsigma = ((1.0/(12.0*_desc[i.idx().p].dphi())))*(-radius(i.idx().p, i.idx().i, i.idx().j+2)
                                                                   + radius(i.idx().p, i.idx().i, i.idx().j-2)
                                                              + 8.0*(radius(i.idx().p, i.idx().i, i.idx().j+1)
                                                                   - radius(i.idx().p, i.idx().i, i.idx().j-1)));
               
               return (1.*sqrt(r*r*pow(cos(rho),-2.)*pow(cos(sigma),-2.)*pow(1.*pow(cos(sigma),-2.)
                      +pow(tan(rho),2.),-4.)*(r*r*(1.*pow(cos(rho),-4.)*pow(cos(sigma),-4.)-pow(cos(rho),-2.)
                      *pow(cos(sigma),-2.)*pow(tan(rho),2.)*pow(tan(sigma),2.))+pow(1.*pow(cos(sigma),-2.)
                      +pow(tan(rho),2.),2.)*(1.*(drdsigma*drdsigma)*pow(cos(rho),-2.)
                      +1.*(drdrho*drdrho)*pow(cos(sigma),-2.)+2.*drdrho*drdsigma*tan(rho)*tan(sigma)))));
            }
            
            
            void output (io_base& io) const;
            
            void input (io_base& io) { }
            
            iter begin() { return iter(this); }
            const_iter begin() const { return const_iter(this); }
            
            /// iterator class to traverse the grid of the slice
            class const_iter
            {
               public :
                        const_iter(spheredata_6patch const * const spheredata_) 
                           : _spheredata(spheredata_),
                             _idx(0,0,0,0,spheredata_->_desc[0].origin_theta() + (spheredata_->lbnd(0)[0]*spheredata_->delta()[0]), 
                                          spheredata_->_desc[0].origin_phi() + (spheredata_->lbnd(0)[1]*spheredata_->delta()[1])) 
                        { }
                        virtual ~const_iter() { }
                        
                        /// dereferencing operator as data-point accessor
                        T operator*() const 
                        { 
                           assert(_spheredata->data.size() == npatches);
                           assert(_spheredata->data[_idx.p].size() > 0);
                           return (*_spheredata).data[_idx.p][_idx.ij]; 
                        }
                        
                        spheredata_6patch<T> const * spheredata() const
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
                              _idx.theta = _spheredata->_desc[0].origin_theta() + (_spheredata->lbnd(_idx.p)[0]*_spheredata->delta()[0]);
                              if (_idx.j >= _spheredata->lsh(_idx.p)[1]-1)
                              {
                                 _idx.p++;
                                 if (_idx.p == _spheredata->npatches)
                                    return;
                                 _idx.theta = _spheredata->_desc[0].origin_theta() + (_spheredata->lbnd(_idx.p)[0]*_spheredata->delta()[0]);
                                 _idx.phi = _spheredata->_desc[0].origin_phi() + (_spheredata->lbnd(_idx.p)[1]*_spheredata->delta()[1]);
                                 _idx.ij = 0;
                                 _idx.i = 0;
                                 _idx.j = 0;
                              }
                              else
                              {
                                 _idx.j++;
                                 _idx.ij++;
                                 _idx.phi+=_spheredata->delta()[1];
                              }
                           }
                           else
                           {
                              _idx.i++;
                              _idx.ij++;
                              _idx.theta+=_spheredata->delta()[0];
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
                           if (_idx.i < _spheredata->nghosts() || 
                               _idx.i > _spheredata->lsh(_idx.p)[0]-_spheredata->nghosts()-1 ||
                               _idx.j < _spheredata->nghosts() || 
                               _idx.j > _spheredata->lsh(_idx.p)[1]-_spheredata->nghosts()-1)
                               return true;
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
                        spheredata_6patch<T> const * const _spheredata;
            };
            
            class iter : public const_iter 
            {
               public :
                        iter(spheredata_6patch* const spheredata_) : _spheredata(spheredata_), const_iter(spheredata_) { }
                        virtual ~iter() { }
                        
                        /// dereferencing operator as data-point accessor
                        T& operator*() 
                        { 
                           assert(_spheredata->data.size() == npatches);
                           assert(_spheredata->data[this->_idx.p].size() > 0);
                           return (*_spheredata).data[this->_idx.p][this->_idx.ij];  
                        }
                        
                        spheredata_6patch<T>* spheredata()
                        {
                           return _spheredata; 
                        }
               private :
                        spheredata_6patch<T>* const _spheredata;
            };
            
            
            /// integrator class
            class integrator
            {
               public :
                        integrator(const spheredata_6patch<T>& sph_) : sph(sph_) { init(); }
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
                           _result *= sph._desc[0].dtheta()*sph._desc[0].dphi();
                           
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
                        void sum(const_iter& it, CCTK_REAL f = 1.0, CCTK_REAL det = 1.0) 
                        {
                           CCTK_REAL weight_theta = 1.0, weight_phi = 1.0;
                        
                           // dont't integrate in ghostzones
                           if (it.ghostzone())
                              return;
                           
                           if (it.idx().i == sph.nghosts())
                           {
                              if (sph.lbnd(it.idx().p)[0] == 0)
                                 weight_theta = 0.5;
                              theta_sum = 0;
                           }
                           if (it.idx().i+sph.lbnd(it.idx().p)[0] == sph.gsh(it.idx().p)[0]-sph.nghosts()-1)
                              weight_theta = 0.5;
                           
                           if (it.idx().j+sph.lbnd(it.idx().p)[1] == sph.nghosts() || it.idx().j+sph.lbnd(it.idx().p)[1] == sph.gsh(it.idx().p)[1]-sph.nghosts()-1)
                              weight_phi = 0.5;
                           
                           theta_sum += *it * det * f
                                     * weight_theta * weight_phi;
                           
                           if (it.idx().i == sph.lsh(it.idx().p)[0]-sph.nghosts()-1)
                              _result += theta_sum;
                        }
               private :
                        const spheredata_6patch<T>& sph;
                        CCTK_REAL theta_sum;
                        CCTK_REAL _result;
                        bool finalized;
            };
            
   private :
   
            /// return decomposed local piece of data for this MPI process
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
               
               for (int p=0; p < npatches; ++p)
               {
                  if (this_proc < used_procs)
                  {
                     _desc[p].lsh() = (_desc[p].gsh() - vect<int,2>(2*spheredata<T>::nghosts(), 2*spheredata<T>::nghosts())) / dims + 2*spheredata<T>::nghosts();
                     
                     vect<int,2> mod = (_desc[p].gsh() - vect<int,2>(2*spheredata<T>::nghosts(), 2*spheredata<T>::nghosts())) % dims;
                  
                     for (int i=0; i < 2; ++i)
                        if (proc_vect[i] < mod[i])
                           _desc[p].lsh()[i]++;
                     
                     assert(_desc[p].lsh()[0] != 2*spheredata<T>::nghosts() && _desc[p].lsh()[1] != 2*spheredata<T>::nghosts());
                     
                     for (int i=0; i < 2; i++)
                     {
                        if (proc_vect[i] == 0)
                           _desc[p].lbnd()[i] = 0;
                        else
                        {
                           _desc[p].lbnd()[i] = _desc[p].lsh()[i] * (proc_vect[i]) - proc_vect[i]*2*spheredata<T>::nghosts();
                           if (proc_vect[i] >= mod[i])
                              _desc[p].lbnd()[i] += mod[i];
                        }
                        
                        if (proc_vect[i] == dims[i]-1)
                           _desc[p].ubnd()[i] = _desc[p].gsh()[i]-1;
                        else
                           _desc[p].ubnd()[i] = _desc[p].lbnd()[i] + _desc[p].lsh()[i]-1;// + spheredata<T>::nghosts();
                     }
         
                     /*for (int i=0; i < 2; i++)
                     {
                        if (proc_vect[i] == dims[i]-1)
                        {
                           _desc[p].lsh()[i] += (_desc[p].gsh()[i] - 2*spheredata<T>::nghosts()) % dims[i];
                           _desc[p].ubnd()[i] += (_desc[p].gsh()[i] - 2*spheredata<T>::nghosts()) % dims[i];
                        }
                     }*/
                     
                     //cout << "----------- " << this_proc << ": " << _lsh << ",  " << _lbnd << ",   " << _ubnd << endl;
                     
                     continue;
                  }
                  else
                  {
                     // this processor carries no data....
                     _desc[p].lsh() = vect<int,2>(0, 0);
                     _desc[p].lbnd() = vect<int,2>(0, 0);
                     _desc[p].ubnd() = vect<int,2>(0, 0);
                  }
               }
            }
            
            /// the data of each patch (comprises the whole sphere)
            vector<vector<T> > data;
            
            /// if we have a non-constant sphere, we will store the
            /// surface radii for each point here
            vector<vector<CCTK_REAL> >* _radii;
            
            /// descriptor for each patch
            class patch_desc
            {
               public :
                        patch_desc() : _gsh(0), _lsh(0), _ubnd(0), _lbnd(0) { }
                        patch_desc(const vect<int, 2>& gsh_,
                                   const vect<int, 2>& lsh_,
                                   const vect<int, 2>& lbnd_,
                                   const vect<int, 2>& ubnd_,
                                   const int nghosts_) 
                           : _nghosts(nghosts_), _gsh(gsh_), _lsh(lsh_), _ubnd(ubnd_), _lbnd(lbnd_),
                             _dtheta((PI/2.0)/(gsh_[0]-1-2*nghosts_)), _dphi((PI/2.0)/(gsh_[1]-1-2*nghosts_))
                        {
                           _origin_theta = -PI/4.0-nghosts_*_dtheta;
                           _origin_phi   = -PI/4.0-nghosts_*_dphi;
                        }
                        
                        virtual ~patch_desc() { }
                        
                        
                        vect<int, 2> gsh() const { return _gsh; }
                        vect<int, 2> lsh() const { return _lsh; }
                        vect<int, 2> lbnd() const { return _lbnd; }
                        vect<int, 2> ubnd() const { return _ubnd; }
                        
                        vect<int, 2>& gsh() { return _gsh; }
                        vect<int, 2>& lsh() { return _lsh; }
                        vect<int, 2>& lbnd() { return _lbnd; }
                        vect<int, 2>& ubnd() { return _ubnd; }
                        
                        CCTK_REAL dtheta() const { return _dtheta; }
                        CCTK_REAL dphi() const { return _dphi; }
                        
                        CCTK_REAL origin_theta() const { return _origin_theta; }
                        CCTK_REAL origin_phi() const { return _origin_phi; }
                        
               private :
                        vect<int, 2> _gsh;
                        vect<int, 2> _lsh;
                        vect<int, 2> _lbnd;
                        vect<int, 2> _ubnd;
                        
                        /// the origin of the global angular coordinates for this patch
                        CCTK_REAL _origin_theta;
                        CCTK_REAL _origin_phi;
                        
                        /// the spacing on the sphere for this patch
                        CCTK_REAL _dphi, _dtheta;
                        
                        /// number of ghostpoints
                        int _nghosts;
            };
            
            vector<patch_desc> _desc;
            
            /// this flag is set initialy according to
            /// whether Llama is activated, the origin is 0, radius=const,
            /// and whether this slice lies on a Llama radial gridpoint.
            bool _uses_Llama;
};


template <typename T>
inline io_base& operator<< (io_base& io, const spheredata_6patch<T>& sd)
{
   sd.output(io);
   
   return io;
} 


}


#endif

