
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

#ifndef _FUNCTIONS_sYlm_
#define _FUNCTIONS_sYlm_


#include <math.h>
#include <cassert>


namespace SPS {

/*
   Calculates the real and imaginary spin-weighted spherical harmonics 
   with given spin s for l,m in spherical coordinates.
   For s=0, we obtain the standard spherical harmonics Ylm.
*/
extern void sYlm(const int ss, const int ll, const int mm, const double theta, const double phi, double* const res_re, double* const res_im);


}

#endif
