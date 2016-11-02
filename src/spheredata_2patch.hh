
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

#ifndef _SPHEREDATA2_
#define _SPHEREDATA2_

#include <vector>

#include "spheredata.hh"

namespace SPS {


using namespace std;


/**
   A stereographic 2-patch system covering the sphere.
   This is probably not what we want for now...it's just
   here for later convenience if someone feels like it...
*/
template <typename T>
class spheredata_2patch : public spheredata<T>
{
   public :
            static int const npatches = 2;
   
            spheredata_2patch() { }
            virtual ~spheredata_2patch() { }
   
   private :
            /// the data of the entire sphere
            vector<vector<T> > data;
};


}


#endif

