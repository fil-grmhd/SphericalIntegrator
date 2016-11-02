
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

#ifndef _COMMSTACK_
#define _COMMSTACK_

#include "mpi.h"
#include "cctk.h"
#include <vector>
#include <cassert>

namespace SPS {


using namespace std;



/**
   This class represents a stack of collective 
   MPI_Allreduce commands. Since MPI_Allreduce is
   expensive, we collect all reductions and do it in one
   single call. 
*/
class commstack
{
   public :
            commstack(MPI_Op op, MPI_Comm comm)
               : collective_buffer(0),
                 _op(op), _comm(comm),
                 was_reduced(false)
            { }
            
            virtual ~commstack() { }
            
            /// puts a value to the collective reduction buffer
            void push(CCTK_REAL val);
            
            /// MPI_Allreduce of the collective buffer
            void reduce();
            
            /// read value from collective buffer. This is checked for reduction.
            CCTK_REAL reduced_val(const int i) const { assert(was_reduced); return collective_buffer[i]; }
            
            /// read value from collective buffer. This may or may not be reduced!
            CCTK_REAL buffer_val(const int i) const { return collective_buffer[i]; }
            
   private :
            
            vector<CCTK_REAL> collective_buffer;
            
            MPI_Op   _op;
            MPI_Comm _comm;
            
            bool was_reduced;
};



}


#endif

