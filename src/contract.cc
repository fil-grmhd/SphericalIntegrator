
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

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "slices.hh"


using namespace SPS;


// the slice ids for each slice and precalculated sYlm
vector<vector<int> > sid_sYlm_re(0);
vector<vector<int> > sid_sYlm_im(0);

vector<int> stored_lmin(0);
vector<int> stored_lmax(0);


extern "C" CCTK_INT SphericalSlice_ContractWithAllPrecalcedsYlm(const CCTK_INT si, const CCTK_INT varno, const CCTK_INT timelevel, const CCTK_INT s, const CCTK_INT lmin, const CCTK_INT lmax, CCTK_COMPLEX* const coeffs);
void precalc_sYlm(const CCTK_INT si, const int s, const int lmin, const int lmax);


extern "C" CCTK_COMPLEX SphericalSlice_ContractWithsYlm(const CCTK_INT varno, const CCTK_INT timelevel, const CCTK_INT s, const CCTK_INT l, const CCTK_INT m)
{
   DECLARE_CCTK_PARAMETERS
   
   assert(timelevel >= 0);
   assert(varno >= 0);
   
   CCTK_COMPLEX result;
   
   if (is_1patch(varno))
   {
      assert(INDEX1P(varno) < slices_1patch.slice().size());
      result = slices_1patch(INDEX1P(varno), timelevel).contract_with_sYlm(s, l, m);
   }
   
   if (is_2patch(varno))
      CCTK_WARN(0, "Uh oh....the idea is good but the world isn't ready yet...");
   
   if (is_6patch(varno))
   {
      assert(INDEX6P(varno) < slices_6patch.slice().size());
      result = slices_6patch(INDEX6P(varno), timelevel).contract_with_sYlm(s, l, m);
   }
   
   // return the surface integral
   return result;
}




extern "C" CCTK_INT SphericalSlice_ContractWithAllsYlm(const CCTK_INT varno, const CCTK_INT timelevel, const CCTK_INT s, const CCTK_INT lmin, const CCTK_INT lmax, CCTK_COMPLEX* const coeffs)
{
   DECLARE_CCTK_PARAMETERS
   
   assert(timelevel >= 0);
   assert(varno >= 0);
   
   int mode_dim = (lmax-lmin+1)*(lmax-lmin+1);
   
   stored_lmin.resize(nslices, 0);
   stored_lmax.resize(nslices, -1);
   
   if (is_1patch(varno))
   {
      assert(INDEX1P(varno) < slices_1patch.slice().size());
      
      const int ID = slices_1patch(INDEX1P(varno), timelevel).ID();
      if (precalc_sYlms && s == 0 && (stored_lmin[ID] > lmin || lmax > stored_lmax[ID]))
         precalc_sYlm(ID, s, lmin, lmax);
         
      if (precalc_sYlms && s == 0 && stored_lmin[ID] <= lmin && lmax <= stored_lmax[ID])
      {
         SphericalSlice_ContractWithAllPrecalcedsYlm(ID, varno, timelevel, s, lmin, lmax, coeffs);
         return 0;
      }
      
      if (!slices_1patch(INDEX1P(varno), timelevel).has_constant_radius() && slices_1patch(INDEX1P(varno), timelevel).nghosts() < 2) 
         CCTK_WARN(0, "ghostzone-width too small!");
                  
      vector<spheredata_1patch<CCTK_REAL>::integrator> intsYlm_re = vector<spheredata_1patch<CCTK_REAL>::integrator>(mode_dim, slices_1patch(INDEX1P(varno), timelevel));
      vector<spheredata_1patch<CCTK_REAL>::integrator> intsYlm_im = vector<spheredata_1patch<CCTK_REAL>::integrator>(mode_dim, slices_1patch(INDEX1P(varno), timelevel));
      
      for (spheredata_1patch<CCTK_REAL>::const_iter i=slices_1patch(INDEX1P(varno), timelevel).begin(); !i.done(); ++i)
      {
         // dont't integrate in ghostzones
         if (i.ghostzone())
            continue;
         
         CCTK_REAL r = slices_1patch(INDEX1P(varno), timelevel).radius(i.idx().p, i.idx().i, i.idx().j);
         CCTK_REAL det = slices_1patch(INDEX1P(varno), timelevel).det(i) / (r*r);
         
         int lm = 0;
         for (int l=lmin; l <= lmax; ++l)
         {
            if (l < abs(s))
            {
               lm+=2*l+1;
               continue;
            }
            
            for (int m=-l; m <= l; ++m)
            {
	       double sYlm_re, sYlm_im;
	       sYlm(s,l,m, i.idx().theta, i.idx().phi, &sYlm_re, &sYlm_im);
               intsYlm_re[lm].sum(i, det, sYlm_re);
               intsYlm_im[lm].sum(i, det, sYlm_im);
               lm++;
            }
         }
      }
      
      // create a collective commstack in order to do only one MPI_Allreduce call instead of 2*l^2
      commstack cs(MPI_SUM, MPI_COMM_WORLD);
      
      int lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               intsYlm_re[lm].finalize(&cs); 
               intsYlm_im[lm].finalize(&cs);
            }
            lm++;
         }
      }
      
      // collective reduction (but only if distrib-method is not constant!!)
      if (slices_1patch(INDEX1P(varno), timelevel).distrib_method() != constant)
         cs.reduce();
      
      lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               coeffs[lm] = CCTK_Cmplx(cs.buffer_val(2*(lm-s*s)), -cs.buffer_val(2*(lm-s*s)+1));
            }
            lm++;
         }
      }
   }
   
   if (is_2patch(varno))
      CCTK_WARN(0, "Uh oh....the idea is good but the world isn't ready yet...");
   
   if (is_6patch(varno))
   {
      assert(INDEX6P(varno) < slices_6patch.slice().size());
      
      const int ID = slices_6patch(INDEX6P(varno), timelevel).ID();
      if (precalc_sYlms && s == 0 && (stored_lmin[ID] > lmin || lmax > stored_lmax[ID]))
         precalc_sYlm(ID, s, lmin, lmax);
         
      if (precalc_sYlms && s == 0 && stored_lmin[ID] >= lmin && lmax <= stored_lmax[ID])
      {
         SphericalSlice_ContractWithAllPrecalcedsYlm(ID, varno, timelevel, s, lmin, lmax, coeffs);
         return 0;
      }
      
      if (!slices_6patch(INDEX6P(varno), timelevel).has_constant_radius() && slices_6patch(INDEX6P(varno), timelevel).nghosts() < 2) 
         CCTK_WARN(0, "ghostzone-width too small!");
                  
      vector<spheredata_6patch<CCTK_REAL>::integrator> intsYlm_re = vector<spheredata_6patch<CCTK_REAL>::integrator>(mode_dim, slices_6patch(INDEX6P(varno), timelevel));
      vector<spheredata_6patch<CCTK_REAL>::integrator> intsYlm_im = vector<spheredata_6patch<CCTK_REAL>::integrator>(mode_dim, slices_6patch(INDEX6P(varno), timelevel));
      
      for (spheredata_6patch<CCTK_REAL>::const_iter i=slices_6patch(INDEX6P(varno), timelevel).begin(); !i.done(); ++i)
      {
         // dont't integrate in ghostzones
         if (i.ghostzone())
            continue;
         
         CCTK_REAL r = slices_6patch(INDEX6P(varno), timelevel).radius(i.idx().p, i.idx().i, i.idx().j);
         CCTK_REAL det = slices_6patch(INDEX6P(varno), timelevel).det(i) / (r*r);
         
         int lm = 0;
         for (int l=lmin; l <= lmax; ++l)
         {
            if (l < abs(s))
            {
               lm+=2*l+1;
               continue;
            }
            
            for (int m=-l; m <= l; ++m)
            {
	       double sYlm_re, sYlm_im;
	       sYlm(s,l,m, slices_6patch(INDEX6P(varno), timelevel).coord_spherical(i.idx().p, i.idx().i, i.idx().j)[0],
                           slices_6patch(INDEX6P(varno), timelevel).coord_spherical(i.idx().p, i.idx().i, i.idx().j)[1],
			   &sYlm_re, &sYlm_im);
               intsYlm_re[lm].sum(i, det, sYlm_re);
               intsYlm_im[lm].sum(i, det, sYlm_im);
               lm++;
            }
         }
      }
            
      commstack cs(MPI_SUM, MPI_COMM_WORLD);
      
      int lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               intsYlm_re[lm].finalize(&cs); 
               intsYlm_im[lm].finalize(&cs);
            }
            lm++;
         }
      }
      
      // collective reduction (but only if distrib-method != constant !!!)
      if (slices_6patch(INDEX6P(varno), timelevel).distrib_method() != constant)
         cs.reduce();
      
      lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               coeffs[lm] = CCTK_Cmplx(cs.buffer_val(2*(lm-s*s)), -cs.buffer_val(2*(lm-s*s)+1));
            }
            lm++;
         }
      }
   }
   
   
   return 0;
}




extern "C" CCTK_INT SphericalSlice_ContractWithAllPrecalcedsYlm(const CCTK_INT si, const CCTK_INT varno, const CCTK_INT timelevel, const CCTK_INT s, const CCTK_INT lmin, const CCTK_INT lmax, CCTK_COMPLEX* const coeffs)
{
   DECLARE_CCTK_PARAMETERS
   
   assert(timelevel >= 0);
   assert(varno >= 0);
   
   int mode_dim = (lmax-lmin+1)*(lmax-lmin+1);
   int off = lmin*lmin;

   if (is_1patch(varno))
   {
      assert(INDEX1P(varno) < slices_1patch.slice().size());
      
      if (!slices_1patch(INDEX1P(varno), 0).has_constant_radius() && slices_1patch(INDEX1P(varno), 0).nghosts() < 2) 
         CCTK_WARN(0, "ghostzone-width too small!");
                  
      vector<integrator_1patch> intsYlm_re = vector<integrator_1patch>(mode_dim, slices_1patch(INDEX1P(varno), 0));
      vector<integrator_1patch> intsYlm_im = vector<integrator_1patch>(mode_dim, slices_1patch(INDEX1P(varno), 0));
      
      for (const_iter_1patch i=slices_1patch(INDEX1P(varno), 0).begin(); !i.done(); ++i)
      {
         // dont't integrate in ghostzones
         if (i.ghostzone())
            continue;
         
         // get local starting index on this processor
         vect<int,2> lbnd = slices_1patch(INDEX1P(varno), 0).lbnd(i.idx().p);
         // get global number of gridpoints for this slice
         vect<int,2> gsh = slices_1patch(INDEX1P(varno), 0).gsh(i.idx().p);
         
         CCTK_REAL r = slices_1patch(INDEX1P(varno), 0).radius(i);
         CCTK_REAL det = slices_1patch(INDEX1P(varno), 0).det(i) / (r*r);
         
         int lm = 0;
         int lm2 = 0;
         for (int l=lmin; l <= lmax; ++l)
         {
            int save_lm2 = lm2;
            lm2 += l;
            for (int m=-l; m < 0; ++m)
            {
               if (l >= abs(s))
               {
                  CCTK_REAL val_sYlm_re = slices_1patch(INDEX1P(sid_sYlm_re[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  CCTK_REAL val_sYlm_im = -slices_1patch(INDEX1P(sid_sYlm_im[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  intsYlm_re[lm].sum(i, det, pow(-1.0, m)*val_sYlm_re);
                  intsYlm_im[lm].sum(i, det, pow(-1.0, m)*val_sYlm_im);
               }
               lm++;
               lm2--;
            }
            lm2 = save_lm2;
            for (int m=0; m <= l; ++m)
            {
               if (l >= abs(s))
               {
                  CCTK_REAL val_sYlm_re = slices_1patch(INDEX1P(sid_sYlm_re[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  CCTK_REAL val_sYlm_im = slices_1patch(INDEX1P(sid_sYlm_im[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  intsYlm_re[lm].sum(i, det, val_sYlm_re);
                  intsYlm_im[lm].sum(i, det, val_sYlm_im);
               }
               lm++;
               lm2++;
            }
         }
      }
      
      // create a collective commstack in order to do only one MPI_Allreduce call instead of 2*l^2
      
      commstack cs(MPI_SUM, MPI_COMM_WORLD);
      
      int lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
               intsYlm_re[lm].finalize(&cs); 
               intsYlm_im[lm].finalize(&cs);
            }
            lm++;
         }
      }
      
      // collective reduction (but only if distrib-method is non-constant!!)
      if (slices_1patch(INDEX1P(varno), timelevel).distrib_method() != constant)
         cs.reduce();
      
      // write back results
      lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
               coeffs[lm] = CCTK_Cmplx(cs.reduced_val(2*(lm-s*s)), -cs.reduced_val(2*(lm-s*s)+1));
            }
            lm++;
         }
      }
      
      
   }
   
   
   if (is_6patch(varno))
   {
      assert(INDEX6P(varno) < slices_6patch.slice().size());
      
      if (!slices_6patch(INDEX6P(varno), 0).has_constant_radius() && slices_6patch(INDEX6P(varno), 0).nghosts() < 2) 
         CCTK_WARN(0, "ghostzone-width too small!");
                  
      vector<integrator_6patch> intsYlm_re = vector<integrator_6patch>(mode_dim, slices_6patch(INDEX6P(varno), 0));
      vector<integrator_6patch> intsYlm_im = vector<integrator_6patch>(mode_dim, slices_6patch(INDEX6P(varno), 0));
      
      for (const_iter_6patch i=slices_6patch(INDEX6P(varno), 0).begin(); !i.done(); ++i)
      {
         // dont't integrate in ghostzones
         if (i.ghostzone())
            continue;
         
         // get local starting index on this processor
         vect<int,2> lbnd = slices_6patch(INDEX6P(varno), 0).lbnd(i.idx().p);
         
         CCTK_REAL r = slices_6patch(INDEX6P(varno), 0).radius(i);
         CCTK_REAL det = slices_6patch(INDEX6P(varno), 0).det(i) / (r*r);
         
         int lm = 0;
         int lm2 = 0;
         for (int l=lmin; l <= lmax; ++l)
         {
            int save_lm2 = lm2;
            lm2 += l;
            for (int m=-l; m < 0; ++m)
            {
               if (l >= abs(s))
               {
                  CCTK_REAL val_sYlm_re = slices_6patch(INDEX6P(sid_sYlm_re[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  CCTK_REAL val_sYlm_im = -slices_6patch(INDEX6P(sid_sYlm_im[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  intsYlm_re[lm].sum(i, det, pow(-1.0, m)*val_sYlm_re);
                  intsYlm_im[lm].sum(i, det, pow(-1.0, m)*val_sYlm_im);
               }
               lm++;
               lm2--;
            }
            lm2 = save_lm2;
            for (int m=0; m <= l; ++m)
            {
               if (l >= abs(s))
               {
                  CCTK_REAL val_sYlm_re = slices_6patch(INDEX6P(sid_sYlm_re[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  CCTK_REAL val_sYlm_im = slices_6patch(INDEX6P(sid_sYlm_im[si][lm2+off]), 0)(i.idx().p, i.idx().i+lbnd[0], i.idx().j+lbnd[1]);
                  intsYlm_re[lm].sum(i, det, val_sYlm_re);
                  intsYlm_im[lm].sum(i, det, val_sYlm_im);
               }
               lm++;
               lm2++;
            }
         }
      }
      
      
      // create a collective commstack in order to do only one MPI_Allreduce call instead of 2*l^2
      
      commstack cs(MPI_SUM, MPI_COMM_WORLD);
      
      int lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
               intsYlm_re[lm].finalize(&cs); 
               intsYlm_im[lm].finalize(&cs);
            }
            lm++;
         }
      }
      
      // collective reduction (but only if distrib-method is non-constant!!!)
      if (slices_6patch(INDEX6P(varno), timelevel).distrib_method() != constant)
         cs.reduce();
      
      // write back results
      lm = 0;
      for (int l=lmin; l <= lmax; ++l)
      {
         for (int m=-l; m <= l; ++m)
         {
            if (l < abs(s))
               coeffs[lm] = CCTK_Cmplx(0.0, 0.0);
            else
            {
               //decomp_vars[lm + mode_dim*(si + nslices*j)] = CCTK_Cmplx(intsYlm_re[lm].finalize(), -intsYlm_im[lm].finalize());
               coeffs[lm] = CCTK_Cmplx(cs.reduced_val(2*(lm-s*s)), -cs.reduced_val(2*(lm-s*s)+1));
            }
            lm++;
         }
      }
   }

   return 0;
}




void precalc_sYlm(const CCTK_INT si, const int s, const int lmin, const int lmax)
{
   DECLARE_CCTK_PARAMETERS

   int mode_dim = (lmax-lmin+1)*(lmax-lmin+1);
   int abs_mode_dim = (lmax+1)*(lmax+1);
   int off = lmin*lmin;
   
   sid_sYlm_re.resize(nslices);
   sid_sYlm_im.resize(nslices);
   
   ostringstream str;
   str << "Precalculating sYlm's for slice " << si;
   CCTK_INFO(str.str().c_str());
   
   for (int si=0; si < nslices; ++si)
   {
      sid_sYlm_re[si].resize(abs_mode_dim, 0);
      sid_sYlm_im[si].resize(abs_mode_dim, 0);
   
      int lm = 0;
      for (int l=lmin; l <= lmax; ++l)
         for (int m=0; m <= l; ++m)  // we only need to store +m modes because -m can be recovered from them via CC!
         {
            // register
            sid_sYlm_re[si][lm+off] = SphericalSlice_RegisterVariable("sYlm_re", si, 1, "const");
            sid_sYlm_im[si][lm+off] = SphericalSlice_RegisterVariable("sYlm_im", si, 1, "const");
            
            // precalculate
            if (is_1patch(sid_sYlm_re[si][lm+off]))
            {
               
               for (iter_1patch i=slices_1patch(INDEX1P(sid_sYlm_re[si][lm+off]), 0).begin(); !i.done(); ++i)
               {
                  double sYlm_re, sYlm_im;
                  sYlm(s,l,m, i.idx().theta, i.idx().phi,
                        &sYlm_re, &sYlm_im);
                  *i = sYlm_re;
                  slices_1patch(INDEX1P(sid_sYlm_im[si][lm+off]), 0)(i) = sYlm_im;
               }
               
            }
            
            if (is_6patch(sid_sYlm_re[si][lm+off]))
            {
               
               for (iter_6patch i=slices_6patch(INDEX6P(sid_sYlm_re[si][lm+off]), 0).begin(); !i.done(); ++i)
               {
                  double sYlm_re, sYlm_im;
                  sYlm(s,l,m, slices_6patch(INDEX6P(sid_sYlm_re[si][lm+off]), 0).coord_spherical(i)[0],
                              slices_6patch(INDEX6P(sid_sYlm_re[si][lm+off]), 0).coord_spherical(i)[1],
                              &sYlm_re, &sYlm_im);
                  *i = sYlm_re;
                  slices_6patch(INDEX6P(sid_sYlm_im[si][lm+off]), 0)(i) = sYlm_im;
               }
               
            }
            
            lm++;
         }
      
   }
   
   stored_lmin[si] = lmin;
   stored_lmax[si] = lmax;
   
}


