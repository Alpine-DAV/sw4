// -*-c++-*-
//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 

#include "Require.h"

#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <fcntl.h>
#include "EW.h"
#include "MaterialSfile.h"
#include "SfileHDF5.h"
#include "Byteswapper.h"

using namespace std;


//-----------------------------------------------------------------------
MaterialSfile::MaterialSfile( EW* a_ew, const string a_file,
			      const string a_directory, bool read_hdf5, bool write_hdf5,
            int horizontalInterval, vector<double> vec_depths):
   mEW(a_ew),
   m_model_file(a_file),
   m_model_dir(a_directory),
   m_read_hdf5(read_hdf5),
   m_write_hdf5(write_hdf5),
   m_horizontalInterval(horizontalInterval),
   m_use_attenuation(a_ew->usingAttenuation())
{
   mCoversAllPoints = false;
   // Check that the depths make sense
   float_sw4 tol = 1e-5;
   float_sw4 depth = 0; // previous depth value, for ascending order
   for (int d=0; d < vec_depths.size(); ++d)
   {
     ASSERT(vec_depths[d] >= (depth+tol));
     depth = vec_depths[d];
     m_vec_depths.push_back(vec_depths[d]);
   }
   if (m_read_hdf5)
     read_sfile();
}

//-----------------------------------------------------------------------
MaterialSfile::~MaterialSfile()
{
}

//-----------------------------------------------------------------------
void MaterialSfile::set_material_properties(std::vector<Sarray> & rhog, 
			       std::vector<Sarray> & csg,
			       std::vector<Sarray> & cpg, 
			       std::vector<Sarray> & xisg, 
			       std::vector<Sarray> & xipg)
{
  if (!m_read_hdf5) // not for reading
    return;
  bool debug=true;
  MPI_Comm comm = MPI_COMM_WORLD;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Barrier(MPI_COMM_WORLD); // FIXME - remove after debugging

  // Assume attenuation arrays defined on all grids if xp defined on grid zero.
  bool use_q = m_use_attenuation && xipg[0].is_defined();
  int ngrids = mEW->mNumberOfGrids;
  int npatch = mMaterial.size();
  for (int g=0; g < ngrids; g++)
  {
    bool curv = mEW->topographyExists() && (g == ngrids-1);
    Sarray& rho=rhog[g];
    Sarray& cp=cpg[g];
    Sarray& cs=csg[g];
    Sarray& qp=xipg[g];
    Sarray& qs=xisg[g];
    if (debug)
    {
      Sarray& data=rho;
      char msg[1000];
      sprintf(msg, "Rank %d, grid %d sw4 bounds [%d:%d,%d:%d,%d:%d]\n",
          myRank, g, data.m_ib, data.m_ie, data.m_jb, data.m_je, 
          data.m_kb, data.m_ke);
      cout << msg;
      cout.flush();
    }

    if (debug)
      cout << "Rank " << myRank << " interpolating grid " << g << endl;
// #pragma omp parallel for
    for (int gj = mEW->m_jStartInt[g]; gj <= mEW->m_jEndInt[g]; ++gj)
    for (int gi = mEW->m_iStartInt[g]; gi <= mEW->m_iEndInt[g]; ++gi)
    {
      float_sw4 gh = mEW->mGridSize[g];
      float_sw4 gx = (gi-1)*gh;
      float_sw4 gy = (gj-1)*gh;
      // float_sw4 gx = min(max(gi-1,0),mEW->m_global_nx[g]*mEW->mGridSize[g];
      // float_sw4 gy = min(max(gj-1,0),mEW->m_global_ny[g])*mEW->mGridSize[g];
      float_sw4 gzmin, gzmax;
      if (curv)
      {
        gzmin = mEW->mZ(gi,gj,mEW->m_kStartInt[g]);
        gzmax = mEW->mZ(gi,gj,mEW->m_kEndInt[g]);
      }
      else
      {
        gzmin = mEW->m_zmin[g] + (mEW->m_kStartInt[g]-1)*mEW->mGridSize[g];
        gzmax = mEW->m_zmin[g] + (mEW->m_kEndInt[g]-1)*mEW->mGridSize[g];
      }

      for (int gk = mEW->m_kStart[g]; gk <= mEW->m_kEnd[g]; ++gk)
      {
        float_sw4 gz;
        if (curv)
          gz = mEW->mZ(gi,gj,gk);
        else
          gz = mEW->m_zmin[g] + (gk-1)*mEW->mGridSize[g];

        bool pasttop = g==(ngrids-1) && gk<1;
        bool pastbot = g==0 && gk>mEW->m_kEndInt[g];
        if (pasttop || pastbot)
          continue; // skip these

        for (int p=0; p < npatch; ++p)
        {
          float_sw4 tol = 1e-2;
          // bottom interface has coarser resolution (except at bottom)
          float ch = (p==0) ? m_hh[0] : m_hh[p-1];
          tol = 1e-2;
          Sarray& ifbot = mInterface[p];
          int cib = (int) floor(1+(gx-m_x0)/ch + tol);
          int cjb = (int) floor(1+(gy-m_y0)/ch + tol);
          cib = min(cib, ifbot.m_ie-1);
          cjb = min(cjb, ifbot.m_je-1);
          float_sw4 zbot00 = -ifbot(cib,cjb,1);
          float_sw4 zbot10 = -ifbot(cib+1,cjb,1);
          float_sw4 zbot01 = -ifbot(cib,cjb+1,1);
          float_sw4 zbot11 = -ifbot(cib+1,cjb+1,1);
          float cwx = (gx-m_x0)/ch - (cib-1); // fraction of ch from cib
          float cwy = (gy-m_y0)/ch - (cjb-1); // fraction of ch from cjb
          float_sw4 zbot = (zbot00*(1-cwx)*(1-cwy) + 
            zbot10*cwx*(1-cwy) + zbot01*(1-cwx)*cwy + 
            zbot11*cwx*cwy);

          float_sw4 h = m_hh[p];
          Sarray& iftop = mInterface[p+1];
          int ib = (int) floor(1+(gx-m_x0)/h + tol);
          int jb = (int) floor(1+(gy-m_y0)/h + tol);
          ib = min(ib, iftop.m_ie-1);
          jb = min(jb, iftop.m_je-1);
          // Interpolate zmin, zmax for this xy location from top interface
          float_sw4 ztop00 = -iftop(ib,jb,1);
          float_sw4 ztop10 = -iftop(ib+1,jb,1);
          float_sw4 ztop01 = -iftop(ib,jb+1,1);
          float_sw4 ztop11 = -iftop(ib+1,jb+1,1);
          float wx = (gx-m_x0)/h - (ib-1); // fraction of h from ib
          float wy = (gy-m_y0)/h - (jb-1); // fraction of h from jb
          float_sw4 ztop = (ztop00*(1-wx)*(1-wy) + 
            ztop10*wx*(1-wy) + ztop01*(1-wx)*wy + 
            ztop11*wx*wy);


          // Do the interpolation if it's in range 
          // OR we're at the top and need to extrapolate?
          bool inz = gz>(ztop-tol) && gz<(zbot+tol);
          bool attop = g==(ngrids-1) && gk==1 && p==(npatch-1);
          if (inz || attop)
          {
            int nk = mMaterial[p].m_ke - mMaterial[p].m_kb + 1;
            float_sw4 hv = (zbot - ztop) / (float_sw4) (nk-1);
            int kb = (int) floor(1+(gz - ztop)/hv+tol);
            kb = max(min(kb, nk-1),1);
            float wz = (gz - ztop)/hv - (kb-1);
            if (gz < ztop)
              wz = 0; // Just copy the value, don't extrapolate
            if (gz > zbot)
              wz = 1; // Just copy the value, don't extrapolate
            float_sw4 rhoi[2], cpi[2], csi[2], qpi[2], qsi[2];

            Sarray& data = mMaterial[p];
            bool outofbounds = ((ib < data.m_ib) || (ib+1 > data.m_ie) || 
                          (jb < data.m_jb) || (jb+1 > data.m_je) || 
                          (kb < data.m_kb) || (kb+1 > data.m_ke));
            if (debug && outofbounds)
            {
              char msg[1000];
              sprintf(msg, "Rank %d, grid %d point [%d,%d,%d], from patch %d, out of bounds [%d,%d,%d] in [%d:%d,%d:%d,%d:%d]\n",
                  myRank, g, gi, gj, gk, p, ib, jb, kb,
                  data.m_ib, data.m_ie, data.m_jb, data.m_je, 
                  data.m_kb, data.m_ke);
              cout << msg;
              cout.flush();
              sprintf(msg, "Rank %d --> values gx=%0.2f, gy=%0.2f, gz=%0.2f, hh=%0.2f\n",
                  myRank, gx,gy,gz,h);
              cout << msg;
              cout.flush();
            }
            if (outofbounds)
              continue; // Can't use this patch
#if 0
#endif // #if 0
            for (int k=0; k < 2; ++k)
            {
              rhoi[k] = mMaterial[p](1,ib,jb,kb+k)*(1-wx)*(1-wy)
              + mMaterial[p](1,ib+1,jb,kb+k)*wx*(1-wy)
              + mMaterial[p](1,ib,jb+1,kb+k)*(1-wx)*wy
              + mMaterial[p](1,ib+1,jb+1,kb+k)*wx*wy;
              cpi[k] = mMaterial[p](2,ib,jb,kb+k)*(1-wx)*(1-wy)
              + mMaterial[p](2,ib+1,jb,kb+k)*wx*(1-wy)
              + mMaterial[p](2,ib,jb+1,kb+k)*(1-wx)*wy
              + mMaterial[p](2,ib+1,jb+1,kb+k)*wx*wy;
              csi[k] = mMaterial[p](3,ib,jb,kb+k)*(1-wx)*(1-wy)
              + mMaterial[p](3,ib+1,jb,kb+k)*wx*(1-wy)
              + mMaterial[p](3,ib,jb+1,kb+k)*(1-wx)*wy
              + mMaterial[p](3,ib+1,jb+1,kb+k)*wx*wy;
              if (use_q)
              {
                qpi[k] = mMaterial[p](4,ib,jb,kb+k)*(1-wx)*(1-wy)
                + mMaterial[p](4,ib+1,jb,kb+k)*wx*(1-wy)
                + mMaterial[p](4,ib,jb+1,kb+k)*(1-wx)*wy
                + mMaterial[p](4,ib+1,jb+1,kb+k)*wx*wy;
                qsi[k] = mMaterial[p](5,ib,jb,kb+k)*(1-wx)*(1-wy)
                + mMaterial[p](5,ib+1,jb,kb+k)*wx*(1-wy)
                + mMaterial[p](5,ib,jb+1,kb+k)*(1-wx)*wy
                + mMaterial[p](5,ib+1,jb+1,kb+k)*wx*wy;
              }
            }
            rho(gi,gj,gk) = rhoi[0]*(1-wz)+rhoi[1]*wz;
            cp(gi,gj,gk) = cpi[0]*(1-wz)+cpi[1]*wz;
            cs(gi,gj,gk) = csi[0]*(1-wz)+csi[1]*wz;
            if (use_q)
            {
              qp(gi,gj,gk) = qpi[0]*(1-wz)+qpi[1]*wz;
              qs(gi,gj,gk) = qsi[0]*(1-wz)+qsi[1]*wz;
            }
            // if (debug && gi==1 && gj==1)
            {
              char msg[1000];
              sprintf(msg, "Grid %d point [%d,%d,%d], gz=%0.2f, wz=%0.2f, from patch %d, [%d,%d,%d]\n",
                  g, gi, gj, gk, gz, wz, p, ib, jb, kb); 
              cout << msg;
              cout.flush();
            }
#if 0
#endif // #if 0
            break; // stop looking in patches, go to next gk
          } // hit on gz in patch
        } // p
      } // k
    } // ij
  } // grids
  mEW->communicate_arrays( rhog );
  mEW->communicate_arrays( csg );
  mEW->communicate_arrays( cpg );
  mEW->material_ic( rhog );
  mEW->material_ic( csg );
  mEW->material_ic( cpg );
  if( use_q )
  {
    mEW->communicate_arrays( xisg );
    mEW->communicate_arrays( xipg );
    mEW->material_ic( xisg );
    mEW->material_ic( xipg );
  }
  MPI_Barrier(MPI_COMM_WORLD); // FIXME - remove after debugging
}

//-----------------------------------------------------------------------
void MaterialSfile::read_sfile()
{
#ifdef USE_HDF5
   string filename = m_model_dir + "/" + m_model_file;
   SfileHDF5::read_sfile_material(filename, *mEW, *this, mMaterial, 
      mInterface);

#else
	 cout << "WARNING: sw4 not compiled with hdf5=yes, " <<
    "--> ignoring MaterialSfile::read_sfile, no-op" << endl;
#endif
}

//-----------------------------------------------------------------------
void MaterialSfile::read_topo(const std::string &file, 
    const std::string &path, EW& ew, Sarray& gridElev,
    float_sw4& lon0, float_sw4& lat0, float_sw4& azim, float_sw4& hh)
{
#ifdef USE_HDF5
   // Don't look in results path
   // string filename = path + "/" + file;
   SfileHDF5::read_sfile_topo(file, ew, gridElev, lon0, lat0, azim, hh);
#else
	 cout << "WARNING: sw4 not compiled with hdf5=yes, " <<
    "--> ignoring MaterialSfile::read_topo, no-op" << endl;
#endif
}

//-----------------------------------------------------------------------
void MaterialSfile::write_sfile(MaterialRfile& rfile)
{
#ifdef USE_HDF5
   SfileHDF5::write_sfile(m_model_file, m_model_dir,
       *mEW, rfile, m_vec_depths, m_horizontalInterval);
#else
	 cout << "WARNING: sw4 not compiled with hdf5=yes, " <<
    "--> ignoring MaterialSfile::write_sfile, no-op" << endl;
#endif
}
