
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

#ifndef _INT_COEFF_
#define _INT_COEFF_


/** 
   This file contains integration coefficients for
   the various methods.
*/


namespace SPS {


using namespace std;


const CCTK_REAL pi = 4.0*atan(1.0);


/*
   Integration coefficients for a 8th
   order scheme (as defined by SBP diagonal norms).
   N+1 = total number of points
   i = current point
*/
inline CCTK_REAL coeff_SBP_8th(int i, int N)
{
      
   assert(N >= 15);

   if ((i == 0) || (i == N))   { return 1498139.0/5080320.0; }
   if ((i == 1) || (i == N-1)) { return 1107307.0/725760.0; }
   if ((i == 2) || (i == N-2)) { return 20761.0/80640.0; }
   if ((i == 3) || (i == N-3)) { return 1304999.0/725760.0; }
   if ((i == 4) || (i == N-4)) { return 299527.0/725760.0; }
   if ((i == 5) || (i == N-5)) { return 103097.0/80640.0; }
   if ((i == 6) || (i == N-6)) { return 670091.0/725760.0; }
   if ((i == 7) || (i == N-7)) { return 5127739.0/5080320.0; }
   if ((i > 7)  || (i < N-7))  { return 1.0; }
   else
      return 0.0;
}


/*
   Integration coefficients for a 6th
   order scheme (as defined by SBP diagonal norms).
   N+1 = total number of points
   i = current point
*/
inline CCTK_REAL coeff_SBP_6th(int i, int N)
{
      
   assert(N >= 11);

   if ((i == 0) || (i == N))   { return 13649.0/43200.0; }
   if ((i == 1) || (i == N-1)) { return 12013.0/8640.0; }
   if ((i == 2) || (i == N-2)) { return 2711.0/4320.0; }
   if ((i == 3) || (i == N-3)) { return 5359.0/4320.0; }
   if ((i == 4) || (i == N-4)) { return 7877.0/8640.0; }
   if ((i == 5) || (i == N-5)) { return 43801.0/43200.0; }
   if ((i > 5)  || (i < N-5))  { return 1.0; }
   else
      return 0.0;
}



/*
   Integration coefficients for Gauss quadrature
   along theta-direction for 1-patch systems.
   This should be saved somewhere because the for-loop makes
   this expensive!
*/
inline CCTK_REAL weight(int i, int N)
{
   CCTK_REAL theta = pi*(i+0.5)/N;
   CCTK_REAL sum1 = 0;
   
   for (int l=0; l <= (N-1)/2.0; ++l)
   {
        sum1 += sin((2.0*l+1.0)*theta)/(2.0*l+1.0);
   }
   
   return sum1;
}


/*
   Integration coefficients for a fourth
   order scheme (variant Simpson's rule).
   N+1 = total number of points
   i = current point
*/
inline CCTK_REAL coeff_variant_Simpson(int i, int N)
{
      
   assert(N >= 7);
      
   if ((i == 0) || (i == N))   { return 17.0/48.0; }
   if ((i == 1) || (i == N-1)) { return 59.0/48.0; }
   if ((i == 2) || (i == N-2)) { return 43.0/48.0; }
   if ((i == 3) || (i == N-3)) { return 49.0/48.0; }
   if ((i > 3)  || (i < N-3))  { return 1.0; }
   else
      return 0;
}



/*
   Integration coefficients for a second
   order midpoint scheme (extended midpoint).
   N+1 = total number of points
   i = current point
*/
inline CCTK_REAL coeff_extended_midpoint(int i, int N)
{
   assert(N >= 1);

   return 1.0;
}




}


#endif

