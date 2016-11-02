
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

#include "spheredata_1patch.hh"

namespace SPS {


template <typename T>
void spheredata_1patch<T>::output (io_base& io) const 
{

   if (data.size() == 0) return;

   const int rank = 2;
   vector<int> dims(rank);
   
   dims[0] = lsh(0)[0];
   dims[1] = lsh(0)[1];
   
   attributes attribs;
   
   // set attributes that shall be attached to the dataset that we are going to write
   // (make it compatible to Carpet's output attributes so that we can use existing
   //  visualization modules)
   attribs << attribute<int>("slice_id", spheredata<T>::ID());
   attribs << attribute<string>("name", spheredata<T>::name());
   attribs << attribute<int>("level", 0);
   
   vector<int> itmp;
   vector<CCTK_REAL> ftmp;
   itmp << _lbnd;
   attribs << attribute<vector<int> >("iorigin", itmp);
   ftmp << vect<CCTK_REAL, 2>(_dtheta, _dphi);
   attribs << attribute<vector<CCTK_REAL> >("delta", ftmp);
   ftmp << vect<CCTK_REAL, 2>(_origin_theta, _origin_phi);
   attribs << attribute<vector<CCTK_REAL> >("origin", ftmp);
   if (spheredata<T>::has_constant_radius())
      attribs << attribute<CCTK_REAL>("radius", radius(0, 0, 0));
   /*else
   {
      // write radius according to overloaded output method
      io.write(attribs,
              (string(spheredata<T>::name())+string("::radius")).c_str(),
              rank,
              dims,
              reinterpret_cast<const double*>(&_radii.front())); 
   }*/
   
  
   
   // write data according to overloaded output method
   io.write(attribs,
            string(spheredata<T>::name()).c_str(),
            rank,
            dims,
            reinterpret_cast<const double*>(&data.front())); 
   
}


  
  // Instantiate template
  template class spheredata_1patch<CCTK_REAL>;
  
} // namespace SPS
