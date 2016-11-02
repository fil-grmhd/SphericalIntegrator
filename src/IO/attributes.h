
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

#ifndef _ATTRIBUTES_
#define _ATTRIBUTES_


#include <vector>
#include <queue>


namespace SPS {


using namespace std;


/**
   A dataset attribute
                        */
template <typename vT, typename iT = string>   // vT = value-type, iT = identifier-type
class attribute
{
   public :
            attribute() : _valid(false) { }
            attribute(iT identifier_, const vT value_) : _identifier(identifier_), _value(value_), _valid(true) { }
            virtual ~attribute() { }
            
            vT    value()       const { return _value; }
            iT    identifier()  const { return _identifier; } 
            bool  valid()       const { return _valid; }
            
   private :
            bool _valid;
            vT   _value;
            iT   _identifier; 
};




/**
   A list of dataset attributes
                                 */
class attributes
{
   public :
            attributes() { }
            virtual ~attributes() { }
            
            void append(const attribute<int>& attrib)
            {
               iattribs.push(attrib);
            }
            
            void append(const attribute<CCTK_REAL>& attrib)
            {
               fattribs.push(attrib);
            }
            
            void append(const attribute<vector<int> >& attrib)
            {
               ivect_attribs.push(attrib);
            }
            
            void append(const attribute<vector<CCTK_REAL > >& attrib)
            {
               fvect_attribs.push(attrib);
            }
            
            void append(const attribute<string>& attrib)
            {
               string_attribs.push(attrib);
            }
            
            void get_next(attribute<int>& attrib)
            {
               if (iattribs.size() == 0)
               {
                  attrib = attribute<int>();
                  return;
               }
               attrib = iattribs.front();
               iattribs.pop();
            }
            
            void get_next(attribute<CCTK_REAL>& attrib)
            {
               if (fattribs.size() == 0)
               {
                  attrib = attribute<CCTK_REAL>();
                  return;
               }
               attrib = fattribs.front();
               fattribs.pop();
            }
            
            void get_next(attribute<vector<int> >& attrib)
            {
               if (ivect_attribs.size() == 0)
               {
                  attrib = attribute<vector<int> >();
                  return;
               }
               attrib = ivect_attribs.front();
               ivect_attribs.pop();
            }
            
            void get_next(attribute<vector<CCTK_REAL> >& attrib)
            {
               if (fvect_attribs.size() == 0)
               {
                  attrib = attribute<vector<CCTK_REAL> >();
                  return;
               }
               attrib = fvect_attribs.front();
               fvect_attribs.pop();
            }
            
            void get_next(attribute<string>& attrib)
            {
               if (string_attribs.size() == 0)
               {
                  attrib = attribute<string>();
                  return;
               }
               attrib = string_attribs.front();
               string_attribs.pop();
            }
            
   private :
            queue<attribute<int,         string> > iattribs;
            queue<attribute<CCTK_REAL ,         string> > fattribs;
            queue<attribute<vector<int>, string> > ivect_attribs;
            queue<attribute<vector<CCTK_REAL >, string> > fvect_attribs;
            queue<attribute<string,      string> > string_attribs;
};



template <typename vT, typename iT>
inline attributes& operator<<(attributes& attribs, const attribute<vT, iT>& attrib)
{
   attribs.append(attrib);
   return attribs;
}


template <typename vT, typename iT>
inline attributes& operator>>(attributes& attribs, attribute<vT, iT>& attrib)
{
   attribs.get_next(attrib);
   return attribs;
}





} // namespace






#endif


