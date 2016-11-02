
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

#include "spheredata_6patch.hh"

namespace SPS {



/// access global angular coordinates (standard theta,phi spherical coordinate system)
template <typename T>
vect<CCTK_REAL, 2> spheredata_6patch<T>::coord_spherical(const int p, const int i, const int j) const 
            { 
               assert(p < npatches);
               
               CCTK_REAL rho = coord(p,i,j)[0];
               CCTK_REAL sigma = coord(p,i,j)[1];
               
               switch (p)
               {
                  case PLUS_X:
                     return vect<CCTK_REAL, 2>(acos( tan(rho) / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ), sigma < 0 ? 2*pi+sigma : sigma );
                  case MINUS_X:
                     return vect<CCTK_REAL, 2>(acos( -tan(rho) / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ), sigma + pi);
                  case PLUS_Y:
                     return vect<CCTK_REAL, 2>(acos( tan(rho) / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ), 
                                               sigma == 0 ? pi/2.0 : sigma < 0 ? atan(1.0/tan(sigma)) + pi : atan(1.0/tan(sigma)));
                  case MINUS_Y:
                     return vect<CCTK_REAL, 2>(acos( -tan(rho) / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ), 
                                               sigma == 0 ? 3.0*pi/2.0 : sigma < 0 ? atan(1.0/tan(sigma)) + 2.0*pi : atan(1.0/tan(sigma)) + pi);
                  case PLUS_Z:
                     return vect<CCTK_REAL, 2>(acos( 1.0 / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ),
                                              (sigma == 0 && rho == 0 ? 0 : (sigma > 0 ? atan(1.0/tan(sigma)*tan(rho)) : atan(1.0/tan(sigma)*tan(rho))+pi  )) );
                  case MINUS_Z:
                     return vect<CCTK_REAL, 2>(acos( 1.0 / sqrt(1.0/pow(cos(sigma),2) + pow(tan(rho),2)) ) + pi, //acos( 1.0 / sqrt(1.0/pow(cos(rho),2) + pow(tan(sigma),2)) ), 
                                              (sigma == 0 && rho == 0 ? 0 : (sigma < 0 ? atan(1.0/tan(sigma)*tan(rho)) : atan(1.0/tan(sigma)*tan(rho))+pi  )) );
               }
               
               return POISON_VAL;
            }



/// returns x-coordinate value of local point i,j on patch p
template <typename T>
CCTK_REAL spheredata_6patch<T>::cart_x(const int p, const int i, const int j) const 
            { 
               assert(p < npatches);
               
               CCTK_REAL x;
               
               switch (p)
               {
                  case PLUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] + radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph);
                     break;
                  }
                  case MINUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] - radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph);
                     break;
                  }
                  case PLUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph) * tan_ph;
                     break;
                  }
                  case MINUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph) * tan_ph;
                     break;
                  }
                  case PLUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu) * tan_nu;
                     break;
                  }
                  case MINUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     x = spheredata<T>::origin()[0] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu) * tan_nu;
                     break;
                  }
               }
               
               return x;
            }
            
/// returns y-coordinate value of local point i,j on patch p
template <typename T>
CCTK_REAL spheredata_6patch<T>::cart_y(const int p, const int i, const int j) const 
            { 
               assert(p < npatches); 
               
               CCTK_REAL y;
               
               switch (p)
               {
                  case PLUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] + radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph) * tan_ph;
                     break;
                  }
                  case MINUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] - radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph) * tan_ph;
                     break;
                  }
                  case PLUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph);
                     break;
                  }
                  case MINUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph);
                     break;
                  }
                  case PLUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu) * tan_mu;
                     break;
                  }
                  case MINUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     y = spheredata<T>::origin()[1] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu) * tan_mu;
                     break;
                  }
               }
               
               return y;
            }
            
/// returns z-coordinate value of local point i,j on patch p
template <typename T>
CCTK_REAL spheredata_6patch<T>::cart_z(const int p, const int i, const int j) const 
            { 
               assert(p < npatches); 
               
               CCTK_REAL z;
               
               switch (p)
               {
                  case PLUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] + radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph) * tan_nu;
                     break;
                  }
                  case MINUS_X:
                  {
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] - radius(p, i, j)/sqrt(1.0 + tan_nu*tan_nu + tan_ph*tan_ph) * tan_nu;
                     break;
                  }
                  case PLUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph) * tan_mu;
                     break;
                  }
                  case MINUS_Y:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_ph = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_ph*tan_ph) * tan_mu;
                     break;
                  }
                  case PLUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] + radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu);
                     break;
                  }
                  case MINUS_Z:
                  {
                     CCTK_REAL const tan_mu = tan(coord(p,i,j)[0]);
                     CCTK_REAL const tan_nu = tan(coord(p,i,j)[1]);
                     z = spheredata<T>::origin()[2] - radius(p, i, j)/sqrt(1.0 + tan_mu*tan_mu + tan_nu*tan_nu);
                     break;
                  }
               }
               
               return z;
            }



template <typename T>
void spheredata_6patch<T>::output (io_base& io) const 
{
   const int rank = 2;
   vector<int> dims(rank);
   
   for (int m=0; m < npatches; ++m)
   {
      
      if (data[m].size() == 0) continue;
      
      dims[0] = lsh(m)[0];
      dims[1] = lsh(m)[1];
      
      attributes attribs;
      
      // set attributes that shall be attached to the dataset that we are going to write
      // (make it compatible to Carpet's output attributes so that we can use existing
      //  visualization modules)
      attribs << attribute<int>("patch", m);
      attribs << attribute<int>("slice_id", spheredata<T>::ID());
      attribs << attribute<string>("name", spheredata<T>::name());
      attribs << attribute<int>("level", 0);
      
      stringstream datasetname;
      datasetname << "patch=" << m << "::" << spheredata<T>::name();
      
      vector<int> itmp;
      vector<CCTK_REAL> ftmp;
      itmp << _desc[m].lbnd();
      attribs << attribute<vector<int> >("iorigin", itmp);
      ftmp << vect<CCTK_REAL, 2>(_desc[m].dtheta(), _desc[m].dphi());
      attribs << attribute<vector<CCTK_REAL> >("delta", ftmp);
      ftmp << vect<CCTK_REAL, 2>(_desc[m].origin_theta(), _desc[m].origin_phi());
      attribs << attribute<vector<CCTK_REAL> >("origin", ftmp);
      if (spheredata<T>::has_constant_radius())
         attribs << attribute<CCTK_REAL>("radius", radius(m, 0, 0));
      /*else
      {
         // write radius according to overloaded output method
         io.write(attribs,
               (datasetname.str()+string("::radius")).c_str(),
               rank,
               dims,
               reinterpret_cast<const double*>(&_radii[m].front())); 
      }*/
      
      // write data according to overloaded output method
      io.write(attribs,
               datasetname.str().c_str(),
               rank,
               dims,
               reinterpret_cast<const double*>(&data[m].front())); 
      
   }
}


  
  // Instantiate template
  template class spheredata_6patch<CCTK_REAL>;
  
} // namespace SPS
