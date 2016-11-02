
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

#include "commstack.hh"


namespace SPS {


using namespace std;



void commstack::push(CCTK_REAL val)
{
   collective_buffer.push_back(val);
   was_reduced = false;
}
            

void commstack::reduce()
{
   vector<CCTK_REAL> dummy(collective_buffer.size(), 0); 
   MPI_Allreduce (&collective_buffer.front(), &dummy.front(), collective_buffer.size(), MPI_DOUBLE, _op, _comm );
   collective_buffer = dummy;
   
   was_reduced = true;
}



}



