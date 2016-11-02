
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

#ifndef _IO_HDF5_SPS_
#define _IO_HDF5_SPS_


#include "io.h"
#include "H5Cpp.h"


namespace SPS {


using namespace H5;


/**
   HDF5 IO class.
                                                               */
class io_hdf5 : public io_base
{
   public :
            io_hdf5(const char* basename_, 
                    const io_mode_t io_mode_ = overwrite, 
                    const int n_active_procs_ = 1, 
                    const int proc_id_ = 0) 
               : file(1), io_base(basename_, ".h5", io_mode_, n_active_procs_, proc_id_) { open(); }
            
            virtual ~io_hdf5() { close(); }
   
            /// return a vector with all registered times 
            virtual vector<CCTK_REAL>  get_times() { return _times; };        
            /// return a vector with all registered iterations
            virtual vector<int> get_cycles() { return _cycles; };       
            /// return the total number of timesteps
            virtual int         get_ntimesteps() { return _cycles.size(); };   
            
            /// open the HDF5 file for the first time and eventually ovewrite it. The File handle
            /// will always be closed after the first opening. This is because we want to re-open/close the
            /// file each time we have performed a read/write process so that if the program crashes during other
            /// than write processes, the HDF5-file will not get corrupted because it was ot properly closed.
            virtual void open()
            {
               // check, if the number of output files matches the number of processes...
               
               // ...if we want to read/overwrite, we don't care about that
               if (_io_mode == overwrite)
               {
                  file = vector<H5File*>(_n_active_procs, NULL);
                  //file.assign(_n_active_procs, NULL);
                  
                  try
                  {
                     Exception::dontPrint();
                     
                     // only write to to the file that belongs to the current process
                     file[_proc_id] = new H5File(filename(_proc_id), H5F_ACC_TRUNC);
                  }
                  // catch failure caused by the H5File operations
                  catch( FileIException error )
                  {
                     error.printError();
                     return;
                  }
               }
               // ...but if we want to append, we have to make sure, we use the same number of files!
               // So first, we try to find out the number of files that have been written
               if (_io_mode == read_only || _io_mode == append)
               {
                  // open master file (which is file 0) and read from attribs
                  // how many files are involved
                  try
                  {
                     Exception::dontPrint();
                  
                     file[0] = new H5File(filename(0), H5F_ACC_RDONLY);
                  }
                  // catch failure caused by the H5File operations
                  catch( FileIException error )
                  {
                     // file is not there, ...so we create a new one
                     // in case we want to append
                     if (_io_mode == append)
                     {
                        try
                        {
                           file = vector<H5File*>(_n_active_procs, NULL);
                           Exception::dontPrint();
                           // create file...
                           file[_proc_id] = new H5File(filename(_proc_id), H5F_ACC_TRUNC);
                           // ...and free handle again.
                           delete file[_proc_id];
                           
                        }
                        catch( FileIException error )
                        {
                           error.printError();
                           return;
                        }
                        _n_io_procs = _n_active_procs;
                     }
                     else
                     {
                        error.printError();
                        return;
                     }
                  }
                  
                  if (file[0]) delete file[0];
                  
                  // try to access "GLOBAL PARAMETERS" group and get number of procs
                  // that wrote the files, etc...
                  get_global_attribs();
               }
               
               if (_io_mode == read_only)
               {
                  // open files for reading 
                  try
                  {
                     Exception::dontPrint();
                     
                     for (int i=0; i < _n_io_procs; i++)
                        file[i] = new H5File(filename(i), H5F_ACC_RDONLY);
                  }
                  // catch failure caused by the H5File operations
                  catch( FileIException error )
                  {
                     error.printError();
                     return;
                  }
               }
               
               if (_io_mode == append)
               {
                  assert(_n_active_procs == _n_io_procs);
                  
                  file = vector<H5File*>(_n_active_procs, NULL);
                  //file.assign(_n_active_procs, NULL);
                  
                  // open files for appending
                  try
                  {
                     Exception::dontPrint();
                     
                     file[_proc_id] = new H5File(filename(_proc_id), H5F_ACC_RDWR);
                  }
                  // catch failure caused by the H5File operations
                  catch( FileIException error )
                  {
                     error.printError();
                     return;
                  }
               }
               
               // afterwards close everything since we will read/write open the file each time we want to access it...
               for (int i=0; i < file.size(); i++)
               {
                  if (file[i]) delete file[i];
                  file[i] = NULL;
               }
            }
            
            virtual void close()
            {
               // delete all opened file handles
               for (int i=0; i < file.size(); i++)
               {
                  if (file[i]) delete file[i];
                  file[i] = NULL;
               }
            }
            
            /// read global IO-attribs that are independent of the hierarchy/field such as the number of
            /// timesteps, number of processes, etc. 
            virtual void get_global_attribs() 
            {
                // open for data-reading
               file[0] = new H5File(filename(0), H5F_ACC_RDONLY);
               
               Group group(file[0]->openGroup("Parameters and Global Attributes"));
               
               Attribute attrib(group.openAttribute("nioprocs"));
               int n_io_procs = 0;
               attrib.read(PredType::NATIVE_INT, &_n_io_procs);
               
               attrib = group.openAttribute("n_timesteps");
               attrib.read(PredType::NATIVE_INT, &_n_timesteps);
               
               try
               {
                  Exception::dontPrint();
                  
                  attrib = group.openAttribute("cycles");
                  DataSpace dspace = attrib.getSpace();
                  _cycles.resize(dspace.getSimpleExtentNpoints());
                  attrib.read(PredType::NATIVE_INT, &_cycles.front());
                  
                  attrib = group.openAttribute("times");
                  _times.resize(_cycles.size());
                  attrib.read(PredType::NATIVE_DOUBLE, &_times.front());
               }
               catch( AttributeIException error )
               {
                  // do nothing, if group was not found
               }
               
               // close again
               delete file[0];
               file[0] = NULL;
            }
            
            /// write global data-attributes into file. If they already exist, they will be overwritten.
            virtual void put_global_attribs() 
            {
               // open for data-appending
               file[_proc_id] = new H5File(filename(_proc_id), H5F_ACC_RDWR);
               
               // try to remove group
               try
               {
                  Exception::dontPrint();
                  file[_proc_id]->unlink("Parameters and Global Attributes");
               }
               catch( FileIException error )
               {
                  // do nothing, if group was not found
               }
               
               // create new group
               Group group(file[_proc_id]->createGroup("Parameters and Global Attributes")); 
               
               Attribute H5attr = group.createAttribute("nioprocs", PredType::NATIVE_INT, DataSpace());
               int val = _n_active_procs;
               H5attr.write(PredType::NATIVE_INT, &val);
               
               H5attr = group.createAttribute("n_timesteps", PredType::NATIVE_INT, DataSpace());
               val = _n_timesteps;
               H5attr.write(PredType::NATIVE_INT, &val);
               
               if (_cycles.size() > 0)
               {
                  hsize_t dim = _cycles.size();
                  DataSpace fspace(1, &dim);
                  H5attr = group.createAttribute("cycles", PredType::NATIVE_INT, fspace);
                  H5attr.write(PredType::NATIVE_INT, &_cycles.front());
                  
                  dim = _times.size();
                  fspace = DataSpace(1, &dim);
                  H5attr = group.createAttribute("times", PredType::NATIVE_DOUBLE, fspace);
                  H5attr.write(PredType::NATIVE_DOUBLE, &_times.front());
               }
               
               // close again
               delete file[_proc_id];
               file[_proc_id] = NULL;
            }
            
            /// set the file-handle to the 
            virtual void set_to_timestep(const int timestep_) { };
            
            /// returns the current timestep in the file handle
            virtual int get_timestep() const { return cur_timestep; }
            
            /// create a new timestep in the file and make it active.
            /// If this timestep already exists, then simply set it to it.
            virtual void create_new_timestep(const int cycle, const CCTK_REAL time) 
            {
               int ts = timestep_already_registered(cycle);
               if (ts >= 0)
               {
                  cur_timestep = ts;
                  cur_cycle = _cycles[ts];
                  cur_time = _times[ts];
               }
               else
               {
                  cur_timestep = _n_timesteps; 
                  _n_timesteps++; 
                  cur_cycle = cycle; 
                  cur_time = time; 
                  _cycles.push_back(cycle);
                  _times.push_back(time);
               }
            };
            
            /// write a CCTK_REAL-dataset to the file
            virtual void write(const attributes& attribs,     // a list of attribs that are attached to the dataset
                               const char* name,                                       // name of the dataset
                               const int  rank,                                   // rank of the dataset
                               const vector<int> dims,                            // extension in each dim
                               const double* buffer)                                  // the buffer that shall be written
            {
               // a copy of the attribs
               attributes a = attribs;
               
               // attach time and cycle attributes
               a << attribute<int>("timestep", cur_cycle);
               a << attribute<double>("time", cur_time);
               
               if (!buffer) return;
               
               hsize_t fdims[dims.size()]; 
               
               // HDF5 writes data in C-order!
               for (int d=0; d < rank; d++)
               {
                  fdims[d] = dims[rank-1-d];
                  /*if (fdims[d] == 0)
                  {
                     return; 
                  }*/
               }
               
               // open for data-appending
               file[_proc_id] = new H5File(filename(_proc_id), H5F_ACC_RDWR);
               
               // Create dataspace for the dataset in the file.
               DataSpace fspace(rank, fdims);
         
               // Create dataset and write it into the file that belongs to the current MPI process.
               stringstream out;
               out << "proc=" << _proc_id << "::cycle=" << cur_cycle << "::" << name;
               DataSet dataset(file[_proc_id]->createDataSet(out.str().c_str(), PredType::NATIVE_DOUBLE, fspace));
               dataset.write(buffer, PredType::NATIVE_DOUBLE);
               
               // attach all attributes
               attribute<int> iattrib;
               a >> iattrib;
               while (iattrib.valid())
               {
                  Attribute H5attr = dataset.createAttribute(iattrib.identifier(), PredType::NATIVE_INT, DataSpace());
                  int val = iattrib.value();
                  H5attr.write(PredType::NATIVE_INT, &val);
                  a >> iattrib;
               }
               attribute<CCTK_REAL> fattrib;
               a >> fattrib;
               while (fattrib.valid())
               {
                  Attribute H5attr = dataset.createAttribute(fattrib.identifier(), PredType::NATIVE_DOUBLE, DataSpace());
                  CCTK_REAL val = fattrib.value();
                  H5attr.write(PredType::NATIVE_DOUBLE, &val);
                  a >> fattrib;
               }
               attribute<string> sattrib;
               a >> sattrib;
               while (sattrib.valid())
               {
                  Attribute H5attr = dataset.createAttribute(sattrib.identifier(), StrType(PredType::C_S1, sattrib.value().size()), DataSpace());
                  const char* val = sattrib.value().c_str();
                  H5attr.write(StrType(PredType::C_S1, sattrib.value().size()), val);
                  a >> sattrib;
               }
               attribute<vector<int> > ivect_attrib;
               a >> ivect_attrib;
               while (ivect_attrib.valid())
               {
                  const hsize_t dim = ivect_attrib.value().size();
                  Attribute H5attr = dataset.createAttribute(ivect_attrib.identifier(), PredType::NATIVE_INT, DataSpace(1, &dim));
                  int val[dim];
                  for (int d=0; d < dim; d++)
                     val[d] = ivect_attrib.value()[d];
                  H5attr.write(PredType::NATIVE_INT, val);
                  a >> ivect_attrib;
               }
               attribute<vector<CCTK_REAL> > fvect_attrib;
               a >> fvect_attrib;
               while (fvect_attrib.valid())
               {
                  const hsize_t dim = fvect_attrib.value().size();
                  Attribute H5attr = dataset.createAttribute(fvect_attrib.identifier(), PredType::NATIVE_DOUBLE, DataSpace(1, &dim));
                  CCTK_REAL val[dim];
                  for (int d=0; d < dim; d++)
                     val[d] = fvect_attrib.value()[d];
                  H5attr.write(PredType::NATIVE_DOUBLE, val);
                  a >> fvect_attrib;
               }
               
               // close again
               delete file[_proc_id];
               file[_proc_id] = NULL;
            }
            
            
            /// write a vector<CCTK_REAL>-dataset to the file
            /// in this case, this is done by creating an individual dataset for each component
            virtual void write(const attributes& attribs,     // a list of attribs that are attached to the dataset
                               const char* name,              // name of the dataset
                               const int  rank,
                               const vector<int> dims,
                               const vector<double>* buffer)     // the buffer that shall be written
            {
               if (!buffer) return;
               
               for (int d=0; d < buffer[0].size(); d++)
               {
                  attributes a = attribs;
                  a << attribute<int>("vector-component", d);
                  ostringstream out;
                  out << name << "::d=" << d;
                  
                  // create new double-buffer
                  int size = 1;
                  for (int i=0; i < rank; i++)
                     size *= dims[i];
                  
                  double* buf = new double[size];
                  
                  // copy vector component to new CCTK_REAL-buffer
                  for (int c=0; c < size; c++)
                     buf[c] = buffer[c][d];
                  
                  // write CCTK_REAL-buffer
                  write(a, out.str().c_str(), rank, dims, buf);
                  delete [] buf;
               }
            }
            
            virtual void read(attributes& attribs,     // a list of CCTK_REAL-attribs that are attached to the dataset
                              const char* name,                                 // name of the dataset
                              CCTK_REAL* buffer)                                  // the buffer where data will be stored
            {
               
            }
            
   private :
            /// all file handles for each processor
            vector<H5File*> file;
};








} // namespace






#endif


