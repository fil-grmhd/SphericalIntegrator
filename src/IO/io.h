
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

#ifndef _IO_
#define _IO_


#include "attributes.h"


namespace SPS {


enum io_mode_t { overwrite, append, read_only };



// conversion of fp to 1d C++ vector
template<typename T>
inline vector<T>& operator<<(vector<T>& a, const double& b)
{
   a.resize(1);
   a[0] = b;

   return a; 
}



/**
   Abstract base class that represents an interface for I/O file formats
                                                                          */
class io_base
{
   public :
            io_base(const char* basename_, const char* extension_,
                    const io_mode_t io_mode_ = overwrite, 
                    const int n_active_procs_ = 1, 
                    const int proc_id_ = 0) 
               : _cycles(0), _times(0), cur_timestep(0), cur_cycle(0), cur_time(0), _n_timesteps(0), _proc_id(proc_id_), _n_active_procs(n_active_procs_), _io_mode(io_mode_), proc(".file_"), extension(extension_), basename(basename_)
            { }
            
            virtual ~io_base() { }
   
            /// return a vector with all registered times 
            virtual vector<double>  get_times() = 0;        
            /// return a vector with all registered iterations
            virtual vector<int> get_cycles() = 0;       
            /// return the total number of timesteps
            virtual int         get_ntimesteps() = 0;   
            
            
            virtual void open() = 0;
            virtual void close() = 0;
            
            /// read global IO-attribs that are independent of the hierarchy/field such as the number of
            /// timesteps, number of processes, etc. 
            virtual void get_global_attribs() = 0;
            
            /// write global data-attributes into file
            virtual void put_global_attribs() = 0;
            
            /// set the file-handle to the 
            virtual void set_to_timestep(const int timestep_) = 0;
            
            /// returns the current timestep in the file handle
            virtual int get_timestep() const { return cur_timestep; }
            
            /// create a new timestep in the file and makes it active
            virtual void create_new_timestep(const int cycle, const double time) = 0;
            
            /// write a CCTK_REAL-dataset to the file
            virtual void write(const int cycle,
                               const double time,
                               const attributes& attribs,     // a list of attribs that are attached to the dataset
                               const char* name,              // name of the dataset
                               const int  rank,               // rank of the dataset
                               const vector<int> dims,        // extents of the dataset
                               const double* buffer)              // the buffer that shall be written
            {
               int ts = timestep_already_registered(cycle);
               if (ts < 0)
                  create_new_timestep(cycle, time);
               else
                  set_to_timestep(ts);
                  
               write(attribs, name, rank, dims, buffer);
            }
            
            /// write a vector<CCTK_REAL>-dataset to the file
            virtual void write(const int cycle,
                               const double time,
                               const attributes& attribs,
                               const char* name,
                               const int  rank,
                               const vector<int> dims,
                               const vector<double>* buffer)
            {
               int ts = timestep_already_registered(cycle);
               if (ts < 0)
                  create_new_timestep(cycle, time);
               else
                  set_to_timestep(ts);
                  
               write(attribs, name, rank, dims, buffer);
            }
            
            /// write a CCTK_REAL-dataset to the current timestep of the file
            virtual void write(const attributes& attribs,
                               const char* name,
                               const int  rank,
                               const vector<int> dims,
                               const double* buffer) = 0;
            
            /// write a vector<CCTK_REAL>-dataset to the current timestep of the file
            virtual void write(const attributes& attribs,
                               const char* name,
                               const int  rank,
                               const vector<int> dims,
                               const vector<double>* buffer) = 0;
            
            virtual void read(attributes& attribs,     // a list of CCTK_REAL-attribs that are attached to the dataset
                              const char* name,        // name of the dataset
                              CCTK_REAL* buffer) = 0;         // the buffer where data will be stored
            
            /// assembles the file-name string
            virtual string filename(const int proc_id_)
            {
               ostringstream out;
               out << basename << proc << proc_id_ << extension;
               return out.str();
            }
            
   protected :
            /// basename of the filename
            const string basename;    
            /// file-name extension
            const string extension;   
            /// proc string
            const string proc;        
            io_mode_t _io_mode;
            
            /// current read/write time step
            int cur_timestep;             
            int cur_cycle;
            CCTK_REAL  cur_time;
            
            /// global I/O attribs
            
            /// this MPI process ID 
            int _proc_id;          
            /// number of files / procs when file was written
            int _n_io_procs;       
            /// number of currently active procs
            int _n_active_procs;   
            /// number of timesteps
            int _n_timesteps;      
            
            vector<int> _cycles;
            vector<CCTK_REAL> _times;
            
            /// check if cycle is already registered
            /// returns -1 if not
            /// return the timestep if found
            int timestep_already_registered(const int cycle)
            {
               if (_n_timesteps > 0)
               {
                  if (cur_cycle != cycle)
                  {
                     vector<int> cycles = get_cycles();
                  
                     for (int i=0; i < cycles.size(); i++)
                     {
                        if (cycles[i] == cycle)
                           return i;
                     }
                     
                     return -1;
                  }
                  else
                     return cur_timestep;
               }
               else
                  return -1;
            }
};








} // namespace






#endif


