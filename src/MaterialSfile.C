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
			      const string a_directory, 
                              int horizontalInterval, vector<double> vec_depths):
   mEW(a_ew),
   m_model_file(a_file),
   m_model_dir(a_directory),
   /* m_read_hdf5(read_hdf5), */
   /* m_write_hdf5(write_hdf5), */
   m_horizontalInterval(horizontalInterval),
   m_use_attenuation(false)
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
   if (a_ew != NULL)
     m_use_attenuation = a_ew->usingAttenuation();
   /* if (m_read_hdf5) */
   /*   read_sfile(); */
     read_sfile2();
}

//-----------------------------------------------------------------------
MaterialSfile::~MaterialSfile()
{
}

/* //----------------------------------------------------------------------- */
/* void MaterialSfile::set_material_properties(std::vector<Sarray> & rhog, */ 
/* 			       std::vector<Sarray> & csg, */
/* 			       std::vector<Sarray> & cpg, */ 
/* 			       std::vector<Sarray> & xisg, */ 
/* 			       std::vector<Sarray> & xipg) */
/* { */
/*   /1* if (!m_read_hdf5) // not for reading *1/ */
/*   /1*   return; *1/ */
/*   bool debug=false; */
/*   MPI_Comm comm = MPI_COMM_WORLD; */
/*   int myRank; */
/*   MPI_Comm_rank(MPI_COMM_WORLD, &myRank); */

/*   // Assume attenuation arrays defined on all grids if xp defined on grid zero. */
/*   bool use_q = m_use_attenuation && xipg[0].is_defined(); */
/*   int ngrids = mEW->mNumberOfGrids; */
/*   int npatch = mMaterial.size(); */
/*   for (int g=0; g < ngrids; g++) */
/*   { */
/*     bool curv = mEW->topographyExists() && (g == ngrids-1); */
/*     Sarray& rho=rhog[g]; */
/*     Sarray& cp=cpg[g]; */
/*     Sarray& cs=csg[g]; */
/*     Sarray& qp=xipg[g]; */
/*     Sarray& qs=xisg[g]; */
/*     if (debug) */
/*     { */
/*       Sarray& data=rho; */
/*       char msg[1000]; */
/*       sprintf(msg, "Rank %d, grid %d sw4 bounds [%d:%d,%d:%d,%d:%d], interior [%d:%d,%d:%d,%d:%d]\n", */
/*           myRank, g, data.m_ib, data.m_ie, data.m_jb, data.m_je, */ 
/*           data.m_kb, data.m_ke, mEW->m_iStartInt[g], mEW->m_iEndInt[g], */ 
/*           mEW->m_jStartInt[g], mEW->m_jEndInt[g], */ 
/*           mEW->m_kStartInt[g], mEW->m_kEndInt[g]); */
/*       cout << msg; */
/*       cout.flush(); */
/*     } */

/*     if (debug) */
/*       cout << "Rank " << myRank << " interpolating grid " << g << endl; */
/* #pragma omp parallel for */
/*     for (int gj = mEW->m_jStartInt[g]; gj <= mEW->m_jEndInt[g]; ++gj) */
/*     for (int gi = mEW->m_iStartInt[g]; gi <= mEW->m_iEndInt[g]; ++gi) */
/*     { */
/*       float_sw4 gh = mEW->mGridSize[g]; */
/*       float_sw4 gx = (gi-1)*gh; */
/*       float_sw4 gy = (gj-1)*gh; */
/*       // float_sw4 gx = min(max(gi-1,0),mEW->m_global_nx[g]*mEW->mGridSize[g]; */
/*       // float_sw4 gy = min(max(gj-1,0),mEW->m_global_ny[g])*mEW->mGridSize[g]; */
/*       float_sw4 gzmin, gzmax; */
/*       if (curv) */
/*       { */
/*         gzmin = mEW->mZ(gi,gj,mEW->m_kStartInt[g]); */
/*         gzmax = mEW->mZ(gi,gj,mEW->m_kEndInt[g]); */
/*       } */
/*       else */
/*       { */
/*         gzmin = mEW->m_zmin[g] + (mEW->m_kStartInt[g]-1)*mEW->mGridSize[g]; */
/*         gzmax = mEW->m_zmin[g] + (mEW->m_kEndInt[g]-1)*mEW->mGridSize[g]; */
/*       } */

/*       for (int gk = mEW->m_kStart[g]; gk <= mEW->m_kEnd[g]; ++gk) */
/*       { */
/*         float_sw4 gz; */
/*         if (curv) */
/*           gz = mEW->mZ(gi,gj,gk); */
/*         else */
/*           gz = mEW->m_zmin[g] + (gk-1)*mEW->mGridSize[g]; */

/*         bool pasttop = g==(ngrids-1) && gk<1; */
/*         bool pastbot = g==0 && gk>mEW->m_kEndInt[g]; */
/*         if (pasttop || pastbot) */
/*           continue; // skip these */

/*         for (int p=0; p < npatch; ++p) */
/*         { */
/*           float_sw4 tol = 1e-2; */
/*           // bottom interface has coarser resolution (except at bottom) */
/*           float ch = (p==0) ? m_hh[0] : m_hh[p-1]; */
/*           tol = 1e-2; */
/*           Sarray& ifbot = mInterface[p]; */
/*           int cib = (int) floor(1+(gx-m_x0)/ch + tol); */
/*           int cjb = (int) floor(1+(gy-m_y0)/ch + tol); */
/*           cib = min(cib, ifbot.m_ie-1); */
/*           cjb = min(cjb, ifbot.m_je-1); */
/*           float_sw4 zbot00 = -ifbot(cib,cjb,1); */
/*           float_sw4 zbot10 = -ifbot(cib+1,cjb,1); */
/*           float_sw4 zbot01 = -ifbot(cib,cjb+1,1); */
/*           float_sw4 zbot11 = -ifbot(cib+1,cjb+1,1); */
/*           float cwx = (gx-m_x0)/ch - (cib-1); // fraction of ch from cib */
/*           float cwy = (gy-m_y0)/ch - (cjb-1); // fraction of ch from cjb */
/*           float_sw4 zbot = (zbot00*(1-cwx)*(1-cwy) + */ 
/*             zbot10*cwx*(1-cwy) + zbot01*(1-cwx)*cwy + */ 
/*             zbot11*cwx*cwy); */

/*           float_sw4 h = m_hh[p]; */
/*           Sarray& iftop = mInterface[p+1]; */
/*           int ib = (int) floor(1+(gx-m_x0)/h + tol); */
/*           int jb = (int) floor(1+(gy-m_y0)/h + tol); */
/*           ib = min(ib, iftop.m_ie-1); */
/*           jb = min(jb, iftop.m_je-1); */
/*           // Interpolate zmin, zmax for this xy location from top interface */
/*           float_sw4 ztop00 = -iftop(ib,jb,1); */
/*           float_sw4 ztop10 = -iftop(ib+1,jb,1); */
/*           float_sw4 ztop01 = -iftop(ib,jb+1,1); */
/*           float_sw4 ztop11 = -iftop(ib+1,jb+1,1); */
/*           float wx = (gx-m_x0)/h - (ib-1); // fraction of h from ib */
/*           float wy = (gy-m_y0)/h - (jb-1); // fraction of h from jb */
/*           float_sw4 ztop = (ztop00*(1-wx)*(1-wy) + */ 
/*             ztop10*wx*(1-wy) + ztop01*(1-wx)*wy + */ 
/*             ztop11*wx*wy); */

/*           // Do the interpolation if it's in range */ 
/*           // OR we're at the top and need to copy? */
/*           // OR we're at the bottom and need to copy? */
/*           bool inz = gz>(ztop-tol) && gz<(zbot+tol); */
/*           bool attop = g==(ngrids-1) && gk>=1 && p==(npatch-1) && gz<ztop+tol; */
/*           bool atbot = g==0 && gk<=mEW->m_kEndInt[g] && p==0 && gz>zbot-tol; */
/*           if (inz || attop || atbot) */
/*           { */
/*             // Make sure we're in the interpolation index range */
/*             ib = max(min(ib, mMaterial[p].m_ie-1),mMaterial[p].m_ib); */
/*             jb = max(min(jb, mMaterial[p].m_je-1),mMaterial[p].m_jb); */
/*             wx = (gx-m_x0)/h - (ib-1); // fraction of h from ib */
/*             wy = (gy-m_y0)/h - (jb-1); // fraction of h from jb */
              
/*             int nk = mMaterial[p].m_ke - mMaterial[p].m_kb + 1; */
/*             float_sw4 hv = (zbot - ztop) / (float_sw4) (nk-1); */
/*             int kb = (int) floor(1+(gz - ztop)/hv+tol); */
/*             kb = max(min(kb, nk-1),1); */
/*             float wz = (gz - ztop)/hv - (kb-1); */
/*             if (gz < ztop) */
/*               wz = 0; // Just copy the value, don't extrapolate */
/*             if (gz > zbot) */
/*               wz = 1; // Just copy the value, don't extrapolate */
/*             if ((wx < -tol || wx > 1+tol) || */ 
/*                 (wy < -tol || wy > 1+tol) || */ 
/*                 (wz < -tol || wz > 1+tol)) */
/*             { */
/*               char msg[1000]; */
/*               /1* sprintf(msg, "Rank %d, grid %d point [%d,%d,%d], from patch %d, [%d:%d,%d:%d,%d:%d], wx=%0.3f, wy=%0.3f, wz=%0.3f!\n", *1/ */
/*               /1*     myRank, g, gi, gj, gk, p, ib, ib+1, jb, jb+1, kb, kb+1, *1/ */
/*               /1*     wx, wy, wz); *1/ */
/*               /1* cout << msg; *1/ */
/*               /1* cout.flush(); *1/ */
/*             } */

/*             float_sw4 rhoi[2], cpi[2], csi[2], qpi[2], qsi[2]; */
/*             Sarray& data = mMaterial[p]; */
/*             bool outofbounds = ((ib < data.m_ib) || (ib+1 > data.m_ie) || */ 
/*                           (jb < data.m_jb) || (jb+1 > data.m_je) || */ 
/*                           (kb < data.m_kb) || (kb+1 > data.m_ke)); */
/*             if (debug && outofbounds) */
/*             { */
/*               char msg[1000]; */
/*               sprintf(msg, "Rank %d, grid %d point [%d,%d,%d], from patch %d, out of bounds [%d:%d,%d:%d,%d:%d] in [%d:%d,%d:%d,%d:%d]\n", */
/*                   myRank, g, gi, gj, gk, p, ib, ib+1, jb, jb+1, kb, kb+1, */
/*                   data.m_ib, data.m_ie, data.m_jb, data.m_je, */ 
/*                   data.m_kb, data.m_ke); */
/*               cout << msg; */
/*               cout.flush(); */
/*               sprintf(msg, "Rank %d --> values gx=%0.2f, gy=%0.2f, gz=%0.2f, hh=%0.2f\n", */
/*                   myRank, gx,gy,gz,h); */
/*               cout << msg; */
/*               cout.flush(); */
/*             } */
/*             if (outofbounds) */
/*               continue; // Can't use this patch */
/*             for (int k=0; k < 2; ++k) */
/*             { */
/*               rhoi[k] = mMaterial[p](1,ib,jb,kb+k)*(1-wx)*(1-wy) */
/*            + mMaterial[p](1,ib+1,jb,kb+k)*wx*(1-wy) */
/*            + mMaterial[p](1,ib,jb+1,kb+k)*(1-wx)*wy */
/*            + mMaterial[p](1,ib+1,jb+1,kb+k)*wx*wy; */
/*               cpi[k] = mMaterial[p](2,ib,jb,kb+k)*(1-wx)*(1-wy) */
/*            + mMaterial[p](2,ib+1,jb,kb+k)*wx*(1-wy) */
/*            + mMaterial[p](2,ib,jb+1,kb+k)*(1-wx)*wy */
/*            + mMaterial[p](2,ib+1,jb+1,kb+k)*wx*wy; */
/*               csi[k] = mMaterial[p](3,ib,jb,kb+k)*(1-wx)*(1-wy) */
/*            + mMaterial[p](3,ib+1,jb,kb+k)*wx*(1-wy) */
/*            + mMaterial[p](3,ib,jb+1,kb+k)*(1-wx)*wy */
/*            + mMaterial[p](3,ib+1,jb+1,kb+k)*wx*wy; */
/*               if (use_q) */
/*               { */
/*                 qpi[k] = mMaterial[p](4,ib,jb,kb+k)*(1-wx)*(1-wy) */
/*              + mMaterial[p](4,ib+1,jb,kb+k)*wx*(1-wy) */
/*              + mMaterial[p](4,ib,jb+1,kb+k)*(1-wx)*wy */
/*              + mMaterial[p](4,ib+1,jb+1,kb+k)*wx*wy; */
/*                 qsi[k] = mMaterial[p](5,ib,jb,kb+k)*(1-wx)*(1-wy) */
/*              + mMaterial[p](5,ib+1,jb,kb+k)*wx*(1-wy) */
/*              + mMaterial[p](5,ib,jb+1,kb+k)*(1-wx)*wy */
/*              + mMaterial[p](5,ib+1,jb+1,kb+k)*wx*wy; */
/*               } */
/*             } */
/*             rho(gi,gj,gk) = rhoi[0]*(1-wz)+rhoi[1]*wz; */
/*             cp(gi,gj,gk) = cpi[0]*(1-wz)+cpi[1]*wz; */
/*             cs(gi,gj,gk) = csi[0]*(1-wz)+csi[1]*wz; */
/*             // Tang */
/*             if (rho(gi, gj, gk) > 2690 || rho(gi, gj, gk) < 1590) */ 
/*                 printf("Rho [%d, %d, %d] = %.1f\n", gi, gj, gk, rho(gi,gj,gk)); */
/*             if (cp(gi, gj, gk) > 5750 || cp(gi, gj, gk) < 700) */ 
/*                 printf("Cp  [%d, %d, %d] = %.1f\n", gi, gj, gk, cp(gi,gj,gk)); */
/*             if (cs(gi, gj, gk) > 5750 || cs(gi, gj, gk) < 80) */ 
/*                 printf("Cs  [%d, %d, %d] = %.1f\n", gi, gj, gk, cs(gi,gj,gk)); */
            
/*             if (use_q) */
/*             { */
/*               qp(gi,gj,gk) = qpi[0]*(1-wz)+qpi[1]*wz; */
/*               qs(gi,gj,gk) = qsi[0]*(1-wz)+qsi[1]*wz; */
/*                 if (qp(gi, gj, gk) > 5750 || qp(gi, gj, gk) < 13) */ 
/*                     printf("Qp  [%d, %d, %d] = %.1f\n", gi, gj, gk, qp(gi,gj,gk)); */
/*                 if (qs(gi, gj, gk) > 5750 || qs(gi, gj, gk) < 26) */ 
/*                     printf("Qs  [%d, %d, %d] = %.1f\n", gi, gj, gk, qs(gi,gj,gk)); */
/*             } */
/*             /1* */
/*             if (debug && gi==1 && gj==1) */
/*             { */
/*               char msg[1000]; */
/*               sprintf(msg, "Grid %d point [%d,%d,%d], gz=%0.2f, wz=%0.2f, from patch %d, [%d,%d,%d]\n", */
/*                   g, gi, gj, gk, gz, wz, p, ib, jb, kb); */ 
/*               cout << msg; */
/*               cout.flush(); */
/*             } */
/*             *1/ */
/*             break; // stop looking in patches, go to next gk */
/*           } // hit on gz in patch */
/*         } // p */
/*       } // k */
/*     } // ij */
/*   } // grids */
/*   mEW->communicate_arrays( rhog ); */
/*   mEW->communicate_arrays( csg ); */
/*   mEW->communicate_arrays( cpg ); */
/*   mEW->material_ic( rhog ); */
/*   mEW->material_ic( csg ); */
/*   mEW->material_ic( cpg ); */
/*   if( use_q ) */
/*   { */
/*     mEW->communicate_arrays( xisg ); */
/*     mEW->communicate_arrays( xipg ); */
/*     mEW->material_ic( xisg ); */
/*     mEW->material_ic( xipg ); */
/*   } */
/* } */

/* //----------------------------------------------------------------------- */
/* // Tang copied and modified from MaterialRfile class */
/* void MaterialSfile::set_material_properties(std::vector<Sarray> & rho, */ 
/*                                              std::vector<Sarray> & cs, */
/*                                              std::vector<Sarray> & cp, */ 
/*                                              std::vector<Sarray> & xis, */ 
/*                                              std::vector<Sarray> & xip ) */
/* { */
/* // Assume attenuation arrays defined on all grids if they are defined on grid zero. */
/*    bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined(); */
/*    size_t outside=0, material=0; */
/*    for( int g=0 ; g < mEW->mNumberOfGrids ; g++ ) { */
/*       bool curvilinear = mEW->topographyExists() && g == mEW->mNumberOfGrids-1; */
/*       float_sw4* rhop=rho[g].c_ptr(); */
/*       float_sw4*  csp=cs[g].c_ptr(); */
/*       float_sw4*  cpp=cp[g].c_ptr(); */
/*       float_sw4* qp, *qs; */
/*       if( use_q  ) { */
/* 	 qs = xis[g].c_ptr(); */
/* 	 qp = xip[g].c_ptr(); */
/*       } */
/*       size_t ni=mEW->m_iEnd[g]-mEW->m_iStart[g]+1; */
/*       size_t nj=mEW->m_jEnd[g]-mEW->m_jStart[g]+1; */
/*       ssize_t ofs = -mEW->m_iStart[g]-ni*(mEW->m_jStart[g])-ni*nj*(mEW->m_kStart[g]); */
/* // Tang comment out for debugging */
/* /1* #pragma omp parallel for reduction(+:material,outside) *1/ */
/*       for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k) */
/* 	 for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j) */
/* 	    for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) { */

/*                 if (g == 1 && j == 49 && i == 273 && k == -1) { */
/*                   int testttt = 1; */
/*                 } */
/* 		float_sw4 x = (i-1)*mEW->mGridSize[g]; */
/* 		float_sw4 y = (j-1)*mEW->mGridSize[g]; */
/* 		float_sw4 z; */
/* 		if( curvilinear ) */
/* 		   z = mEW->mZ(i,j,k); */
/* 		else */
/* 		   z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g]; */
/*                 // Tang: (x, y, z) is the coordinate of current grid point */
/* 		size_t ind = ofs + i + ni*j + ni*nj*k; */
/* 		if( inside( x, y, z ) ) { */

/* 		   material++; */
/*                    int gr = m_npatches-1; */
/*                    // Tang: gr is the patch id that has the current sw4 grid point's data */
/*                    // Tang: need to use the interface value to determine which patch the current point is in */
/*                    int i0, j0, k0; */
/* 		   while( gr > 0 ) { */
/*                      i0 = static_cast<int>( trunc( 1 + (x-m_x0)/m_hh[gr] ) ); */
/*                      j0 = static_cast<int>( trunc( 1 + (y-m_y0)/m_hh[gr] ) ); */

/*                      if (z > mInterface[gr](i0*2-1, j0*2-1, 1)) { */
/*                        break; */
/*                      } */
/* 		     gr--; */
/*                    } */

/*                    // debug */
/*                    float_sw4 t1 = mInterface[gr+1](i0, j0, 1); */
/*                    float_sw4 t2 = mInterface[gr](i0*2-1, j0*2-1, 1); */
/*                    float_sw4 t3 = mInterface[gr](i0, j0, 1); */

/*                    float_sw4 tmph; */
/*                    // Top level interface have same grid size */
/*                    // The bottom interface of a sfile patch is 2x smaller than its top interface */
/*                    if (gr == 0) */ 
/*                        tmph = mInterface[1](i0, j0, 1) - mInterface[0](i0, j0, 1); */
/*                    else */ 
/*                        tmph = mInterface[gr+1](i0, j0, 1) - mInterface[gr](i0*2-1, j0*2-1, 1); */

/*                    // Update the current vertical grid height and z-base with the sfile curvilinear grid */
/*                    m_hv[gr] = tmph / (m_nk[gr]-1); */
/*                    m_z0[gr] = mInterface[gr](i0*2-1, j0*2-1, 1); */


/*                    // Tang: we are using curvilinear grid in sfile */ 
/*                    k0 = static_cast<int>( trunc( 1 + (z-m_z0[gr])/m_hv[gr]) ); */

/* 		   // Use bilinear interpolation always: */
/*                    bool intp_cubic = false; */

/* 	   // Bias stencil near the boundary, need to communicate arrays afterwards. */
/* 		   if( i0 <= m_ifirst[gr] ) { */
/* 		      i0 = m_ifirst[gr]; */
/*                       intp_cubic = false; */
/* 		   } */
/* 		   if( i0 >= m_ilast[gr]-1 ) { */
/* 		      i0 = m_ilast[gr]-1; */
/*                       intp_cubic = false; */
/* 		   } */
/* 		   if( j0 <= m_jfirst[gr] ) { */
/* 		      j0 = m_jfirst[gr]; */
/* 		      intp_cubic = false; */
/* 		   } */
/* 		   if( j0 >= m_jlast[gr]-1 ) { */
/* 		      j0 = m_jlast[gr]-1; */
/*                       intp_cubic = false; */
/* 		   } */
/* 		   if( k0 <= m_kfirst[gr] ) { */
/* 		      k0 = m_kfirst[gr]; */
/*                       intp_cubic = false; */
/* 		   } */
/* 		   if( k0 >= m_klast[gr]-1 ) { */
/* 		      k0 = m_klast[gr]-1; */
/*                       intp_cubic = false; */
/* 		   } */

/* 		   // cubic Hermite intp. */
/*                    if( intp_cubic ) { */
/*                        float_sw4 r = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr]; */
/*                        float_sw4 wghx[4] = {0.5*r*(-r*r+2*r-1),0.5*(3*r*r*r-5*r*r+2),0.5*r*(-3*r*r+4*r+1),0.5*r*r*(r-1)}; */
/*                        float_sw4 s = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr]; */
/*                        float_sw4 wghy[4] = {0.5*s*(-s*s+2*s-1),0.5*(3*s*s*s-5*s*s+2),0.5*s*(-3*s*s+4*s+1),0.5*s*s*(s-1)}; */
/*                        float_sw4 t=  (z-( (k0-1)*m_hv[gr]+m_z0[gr]) )/m_hv[gr]; */
/*                        float_sw4 wghz[4] = {0.5*t*(-t*t+2*t-1),0.5*(3*t*t*t-5*t*t+2),0.5*t*(-3*t*t+4*t+1),0.5*t*t*(t-1)}; */
/*                        rhop[ind]=cpp[ind]=csp[ind]=0; */
/*                        if( use_q ) */
/*                           qp[ind]=qs[ind]=0; */
/*                        for( int kk=0 ; kk < 4 ; kk++ ) */
/*                           for( int jj=0 ; jj < 4 ; jj++ ) */
/*                              for( int ii=0 ; ii < 4 ; ii++ ) */
/*                              { */
/*                                 float_sw4 wgh = wghx[ii]*wghy[jj]*wghz[kk]; */
/*                                 rhop[ind] += wgh*mMaterial_rho[gr](1,i0-1+ii,j0-1+jj,k0-1+kk); */
/*                                 cpp[ind] += wgh*mMaterial_cp[gr](1,i0-1+ii,j0-1+jj,k0-1+kk); */
/*                                 csp[ind] += wgh*mMaterial_cs[gr](1,i0-1+ii,j0-1+jj,k0-1+kk); */
/*                                 if( use_q ) { */
/*                                    qp[ind] += wgh*mMaterial_qp[gr](1,i0-1+ii,j0-1+jj,k0-1+kk); */
/*                                    qs[ind] += wgh*mMaterial_qs[gr](1,i0-1+ii,j0-1+jj,k0-1+kk); */
/*                                 } */
/*                              } */
/* 		   } */
/* 		   else { */
/*    		   // bilinear intp. */
/*                       float_sw4 wghx = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr]; */
/*                       float_sw4 wghy = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr]; */
/*                       float_sw4 wghz = (z-( (k0-1)*m_hv[gr]+m_z0[gr]) )/m_hv[gr]; */

/*                       // Debug */
/*                       if (wghx > 1) */ 
/*                           wghx = 1; */
/*                       if (wghx < 0) */ 
/*                           wghx = 0; */

/*                       if (wghy > 1) */ 
/*                           wghy = 1; */
/*                       if (wghy < 0) */ 
/*                           wghy = 0; */

/*                       if (wghz > 1) */ 
/*                           wghz = 1; */
/*                       if (wghz < 0) */ 
/*                           wghz = 0; */

/*                       rhop[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0) + */ 
/*                                    wghx*mMaterial_rho[gr](1,i0+1,j0,k0) ) + */
/*    			    	   wghy*(     (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0) + */
/*                                    wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0) ) ) + */ 
/*    		                   wghz*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0+1) + */
/*                                    wghx*mMaterial_rho[gr](1,i0+1,j0,k0+1) ) + */
/*    				   wghy*(    (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0+1) + */ 
/*                                    wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0+1) ) ); */

/*                        cpp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0) + */ 
/*                                    wghx*mMaterial_cp[gr](1,i0+1,j0,k0) ) + */
/*                                    wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0) + */ 
/*                                    wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0) ) ) + */ 
/*                                    wghz*(  (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0+1) + */ 
/*                                    wghx*mMaterial_cp[gr](1,i0+1,j0,k0+1) ) + */
/*                                    wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0+1)+ */ 
/*                                    wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0+1) ) ); */
       
/*                        csp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0) + */ 
/*                                    wghx*mMaterial_cs[gr](1,i0+1,j0,k0) ) + */
/*                                    wghy*(    (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0) + */ 
/*                                    wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0) ) ) + */ 
/*                                    wghz*(   (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0+1) + */ 
/*                                    wghx*mMaterial_cs[gr](1,i0+1,j0,k0+1) ) + */
/*                                    wghy*(   (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0+1)+ */ 
/*                                    wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0+1) ) ); */

/*                        if( use_q ) { */
/*                           qp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0) + */ 
/*                                      wghx*mMaterial_qp[gr](1,i0+1,j0,k0) ) + */
/*                                      wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0) + */ 
/*                                      wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0) ) ) + */ 
/*                                      wghz*( (1-wghy)*(   (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0+1) + */ 
/*                                      wghx*mMaterial_qp[gr](1,i0+1,j0,k0+1) ) + */
/*                                      wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0+1)+ */ 
/*                                      wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0+1) ) ); */

/*                           qs[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0) + */ 
/*                                      wghx*mMaterial_qs[gr](1,i0+1,j0,k0) ) + */
/*                                      wghy*(      (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0) + */ 
/*                                      wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0) ) ) + */ 
/*                                      wghz*(   (1-wghy)*(   (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0+1) + */ 
/*                                      wghx*mMaterial_qs[gr](1,i0+1,j0,k0+1) ) + */
/*                                      wghy*(     (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0+1)+ */ 
/*                                      wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0+1) ) ); */
/*                        } */
/*                        float_sw4 rho_tmp = rho[g](i, j, k); */
/*                        // Tang debug */
/*                        /1* ASSERT(rhop[ind] >= 1000); *1/ */
/*                        /1* ASSERT(cpp[ind] >= 600); *1/ */
/*                        /1* ASSERT(csp[ind] >= 70); *1/ */
/*                        /1* if (rhop[ind] > 3000 || rhop[ind] < 1500) *1/ */ 
/*                        /1*     printf("Error with rho[%d]: %.2f\n", ind, rhop[ind]); *1/ */
/*                        /1* if (cpp[ind] > 5700  || cpp[ind] < 600) *1/ */ 
/*                        /1*     printf("Error with cp[%d]: %.2f\n", ind, cpp[ind]); *1/ */
/*                        /1* if (csp[ind] > 3400|| csp[ind] < 70) *1/ */ 
/*                        /1*     printf("Error with cs[%d]: %.2f\n", ind, csp[ind]); *1/ */
/*                        /1* if( use_q ) { *1/ */
/*                        /1*   if (qp[ind] > 377|| qp[ind] < 12) *1/ */ 
/*                        /1*       printf("Error with qp[%d]: %.2f\n", ind, qp[ind]); *1/ */
/*                        /1*   if (qs[ind] > 752|| qs[ind] < 25) *1/ */ 
/*                        /1*       printf("Error with qs[%d]: %.2f\n", ind, qs[ind]); *1/ */
/*                        /1* } *1/ */
     
/*                    } // End else */
/* 		} // End if inside */
/* 		else */
/* 		   outside++; */
/* 	    } */

/*     for (int gk = mEW->m_kStartInt[g]; gk <= mEW->m_kEndInt[g]; ++gk) */
/*     /1* for (int gk = mEW->m_kStart[g]; gk <= mEW->m_kEnd[g]; ++gk) *1/ */
/*       for (int gj = mEW->m_jStartInt[g]; gj <= mEW->m_jEndInt[g]; ++gj) */
/*         for (int gi = mEW->m_iStartInt[g]; gi <= mEW->m_iEndInt[g]; ++gi) { */
/*           if (rho[g](gi,gj,gk) < 0) { */
/*             printf("Error with grid %d rho(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, rho[g](gi,gj,gk), rho[g](gi,gj,gk+1)); */
/*             rho[g](gi,gj,gk) = rho[g](gi,gj,gk+1); */
/*           } */
/*           if (cp[g](gi,gj,gk) < 100) { */
/*             printf("Error with grid %d cp(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, cp[g](gi,gj,gk), cp[g](gi,gj,gk+1)); */
/*             cp[g](gi,gj,gk) = cp[g](gi,gj,gk+1); */
/*           } */
/*           if (cs[g](gi,gj,gk) < 0) { */
/*             printf("Error with grid %d cs(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, cs[g](gi,gj,gk), cs[g](gi,gj,gk+1)); */
/*             cs[g](gi,gj,gk) = cs[g](gi,gj,gk+1); */
/*           } */
/*           if( use_q ) { */
/*             if (xip[g](gi,gj,gk) < 0) { */
/*               printf("Error with grid %d xip(%d, %d, %d): %.2f, changing to %.2f\n", gi, g, gj, gk, xip[g](gi,gj,gk), xip[g](gi,gj,gk+1)); */
/*               xip[g](gi,gj,gk) = xip[g](gi,gj,gk+1); */
/*             } */
/*             if (xis[g](gi,gj,gk) < 0) { */
/*               printf("Error with grid %d xis(%d, %d, %d): %.2f, changing to %.2f\n", gi, g, gj, gk, xis[g](gi,gj,gk), xis[g](gi,gj,gk+1)); */
/*               xis[g](gi,gj,gk) = xis[g](gi,gj,gk+1); */
/*             } */
/*           } */

/*         } */

/*     if (mEW->getRank() == 0) */ 
/*       printf("After interpolation from sfile\n"); */

/*     printf("Rank %d grid %d: rho min = %.2f, max = %.2f\n", mEW->getRank(), g, rho[g].minimum(), rho[g].maximum()); */
/*     printf("Rank %d grid %d: cp min = %.2f, max = %.2f\n", mEW->getRank(), g, cp[g].minimum(), cp[g].maximum()); */
/*     printf("Rank %d grid %d: cs min = %.2f, max = %.2f\n", mEW->getRank(), g, cs[g].minimum(), cs[g].maximum()); */

/*     if( use_q ) { */
/*       printf("Rank %d grid %d: xip min = %.2f, max = %.2f\n", mEW->getRank(), g, xip[g].minimum(), xip[g].maximum()); */
/*       printf("Rank %d grid %d: xis min = %.2f, max = %.2f\n", mEW->getRank(), g, xis[g].minimum(), xis[g].maximum()); */
/*     } */

/*    } // end for g... */

/*    mEW->communicate_arrays( rho ); */
/*    mEW->communicate_arrays( cs ); */
/*    mEW->communicate_arrays( cp ); */
/*    mEW->material_ic( rho ); */
/*    mEW->material_ic( cs ); */
/*    mEW->material_ic( cp ); */
/*    if( use_q ) { */
/*       mEW->communicate_arrays( xis ); */
/*       mEW->communicate_arrays( xip ); */
/*       mEW->material_ic( xis); */
/*       mEW->material_ic( xip ); */
/*    } */

/*    size_t materialSum, outsideSum; */
/*    int mpisizelong, mpisizelonglong, mpisizeint; */
/*    MPI_Type_size(MPI_LONG,&mpisizelong ); */
/*    MPI_Type_size(MPI_LONG_LONG,&mpisizelonglong ); */
/*    MPI_Type_size(MPI_INT,&mpisizeint ); */
/*    if( sizeof(size_t) == mpisizelong ) { */
/*       MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*       MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*    } */
/*    else if( sizeof(size_t) == mpisizelonglong ) { */
/*       MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*       MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*    } */
/*    else if( sizeof(size_t) == mpisizeint ) { */
/*       MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*       MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*    } */
/*    else { */
/*       int materialsumi, outsidesumi, materiali=material, outsidei=outside; */
/*       MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*       MPI_Reduce(&outsidei,   &outsidesumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ); */
/*       materialSum=materialsumi; */
/*       outsideSum=outsidesumi; */
/*    } */
/*    if (mEW->getRank() == 0) */
/*       //      cout << endl */ 
/*       //           << "--------------------------------------------------------------\n" */
/*       //           << "Sfile Initialized Node Types: " << endl */
/*       //           << "   Material:        " << materialSum << endl */
/*       //           << endl */
/*       //           << "*Outside Domain:    " << outsideSum << endl */
/*       //           << endl */ 
/*       //           << "--------------------------------------------------------------\n" */
/*       //           << endl; */
/*       cout << endl */
/* 	   << "sfile command: outside = " << outsideSum << ", material = " << materialSum << endl; */
/* } */


//-----------------------------------------------------------------------
// Tang copied and modified from MaterialRfile class
void MaterialSfile::set_material_properties(std::vector<Sarray> & rho, 
                                             std::vector<Sarray> & cs,
                                             std::vector<Sarray> & cp, 
                                             std::vector<Sarray> & xis, 
                                             std::vector<Sarray> & xip )
{
// Assume attenuation arrays defined on all grids if they are defined on grid zero.
   bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined();
   size_t outside=0, material=0;
   for( int g=0 ; g < mEW->mNumberOfGrids ; g++ ) {
      bool curvilinear = mEW->topographyExists() && g == mEW->mNumberOfGrids-1;
      /* float_sw4* rhop=rho[g].c_ptr(); */
      /* float_sw4*  csp=cs[g].c_ptr(); */
      /* float_sw4*  cpp=cp[g].c_ptr(); */
      /* float_sw4* qp, *qs; */
      /* if( use_q  ) { */
	 /* qs = xis[g].c_ptr(); */
	 /* qp = xip[g].c_ptr(); */
      /* } */
      size_t ni=mEW->m_iEnd[g]-mEW->m_iStart[g]+1;
      size_t nj=mEW->m_jEnd[g]-mEW->m_jStart[g]+1;
// Tang comment out for debugging
/* #pragma omp parallel for reduction(+:material,outside) */
      for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k)
	 for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j)
	    for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i) {

		float_sw4 x = (i-1)*mEW->mGridSize[g];
		float_sw4 y = (j-1)*mEW->mGridSize[g];
		float_sw4 z;
		if( curvilinear )
		   z = mEW->mZ(i,j,k);
		else
		   z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];

                if (g == 1 && ((j == 49 && i == 273 && k == 1) || (j == 48 && i == 274 && k == 1) ||
                               (j == 49 && i == 274 && k == 1) || (j == 49 && i == 275 && k == 1))) {
                    printf("(%d, %d, %d): x=%.2f, y=%.2f, z=%.2f \n", i, j, k, x, y, z);
                }
                // Tang: (x, y, z) is the coordinate of current grid point
		if( inside( x, y, z ) ) {

		   material++;
                   int gr = m_npatches-1;
                   // Tang: gr is the patch id that has the current sw4 grid point's data
                   // Tang: need to use the interface value to determine which patch the current point is in
                   int i0, j0, k0;
		   while( gr > 0 ) {
                     i0 = static_cast<int>( trunc( 1 + (x-m_x0)/m_hh[gr] ) );
                     j0 = static_cast<int>( trunc( 1 + (y-m_y0)/m_hh[gr] ) );

                     if (z > mInterface[gr](i0*2-1, j0*2-1, 1)) {
                       break;
                     }
		     gr--;
                   }

                   // debug
                   float_sw4 t1 = mInterface[gr+1](i0, j0, 1);
                   float_sw4 t2 = mInterface[gr](i0*2-1, j0*2-1, 1);
                   float_sw4 t3 = mInterface[gr](i0, j0, 1);

                   float_sw4 tmph;
                   // Top level interface have same grid size
                   // The bottom interface of a sfile patch is 2x smaller than its top interface
                   if (gr == 0) 
                       tmph = mInterface[1](i0, j0, 1) - mInterface[0](i0, j0, 1);
                   else 
                       tmph = mInterface[gr+1](i0, j0, 1) - mInterface[gr](i0*2-1, j0*2-1, 1);

                   // Update the current vertical grid height and z-base with the sfile curvilinear grid
                   m_hv[gr] = tmph / (m_nk[gr]-1);
                   m_z0[gr] = mInterface[gr](i0*2-1, j0*2-1, 1);


                   // Tang: we are using curvilinear grid in sfile 
                   k0 = static_cast<int>( trunc( 1 + (z-m_z0[gr])/m_hv[gr]) );

		   // Use bilinear interpolation always:
                   bool intp_cubic = false;

	   // Bias stencil near the boundary, need to communicate arrays afterwards.
		   if( i0 <= m_ifirst[gr] ) {
		      i0 = m_ifirst[gr];
                      intp_cubic = false;
		   }
		   if( i0 >= m_ilast[gr]-1 ) {
		      i0 = m_ilast[gr]-1;
                      intp_cubic = false;
		   }
		   if( j0 <= m_jfirst[gr] ) {
		      j0 = m_jfirst[gr];
		      intp_cubic = false;
		   }
		   if( j0 >= m_jlast[gr]-1 ) {
		      j0 = m_jlast[gr]-1;
                      intp_cubic = false;
		   }
		   if( k0 <= m_kfirst[gr] ) {
		      k0 = m_kfirst[gr];
                      intp_cubic = false;
		   }
		   if( k0 >= m_klast[gr]-1 ) {
		      k0 = m_klast[gr]-1;
                      intp_cubic = false;
		   }

		   // cubic Hermite intp.
                   if( intp_cubic ) {
                       float_sw4 r = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr];
                       float_sw4 wghx[4] = {0.5*r*(-r*r+2*r-1),0.5*(3*r*r*r-5*r*r+2),0.5*r*(-3*r*r+4*r+1),0.5*r*r*(r-1)};
                       float_sw4 s = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr];
                       float_sw4 wghy[4] = {0.5*s*(-s*s+2*s-1),0.5*(3*s*s*s-5*s*s+2),0.5*s*(-3*s*s+4*s+1),0.5*s*s*(s-1)};
                       float_sw4 t=  (z-( (k0-1)*m_hv[gr]+m_z0[gr]) )/m_hv[gr];
                       float_sw4 wghz[4] = {0.5*t*(-t*t+2*t-1),0.5*(3*t*t*t-5*t*t+2),0.5*t*(-3*t*t+4*t+1),0.5*t*t*(t-1)};
                       rho[g](i, j, k)=cp[g](i, j, k)=cs[g](i, j, k)=0;
                       if( use_q )
                          xip[g](i, j, k)=xis[g](i, j, k)=0;
                       for( int kk=0 ; kk < 4 ; kk++ )
                          for( int jj=0 ; jj < 4 ; jj++ )
                             for( int ii=0 ; ii < 4 ; ii++ )
                             {
                                float_sw4 wgh = wghx[ii]*wghy[jj]*wghz[kk];
                                rho[g](i, j, k) += wgh*mMaterial_rho[gr](1,i0-1+ii,j0-1+jj,k0-1+kk);
                                cp[g](i, j, k) += wgh*mMaterial_cp[gr](1,i0-1+ii,j0-1+jj,k0-1+kk);
                                cs[g](i, j, k) += wgh*mMaterial_cs[gr](1,i0-1+ii,j0-1+jj,k0-1+kk);
                                if( use_q ) {
                                   xip[g](i, j, k) += wgh*mMaterial_qp[gr](1,i0-1+ii,j0-1+jj,k0-1+kk);
                                   xis[g](i, j, k) += wgh*mMaterial_qs[gr](1,i0-1+ii,j0-1+jj,k0-1+kk);
                                }
                             }
		   }
		   else {
   		   // bilinear intp.
                      float_sw4 wghx = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr];
                      float_sw4 wghy = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr];
                      float_sw4 wghz = (z-( (k0-1)*m_hv[gr]+m_z0[gr]) )/m_hv[gr];

                      // Debug
                      if (wghx > 1) 
                          wghx = 1;
                      if (wghx < 0) 
                          wghx = 0;

                      if (wghy > 1) 
                          wghy = 1;
                      if (wghy < 0) 
                          wghy = 0;

                      if (wghz > 1) 
                          wghz = 1;
                      if (wghz < 0) 
                          wghz = 0;

                      rho[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0) + 
                                   wghx*mMaterial_rho[gr](1,i0+1,j0,k0) ) +
   			    	   wghy*(     (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0) +
                                   wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0) ) ) + 
   		                   wghz*( (1-wghy)*( (1-wghx)*mMaterial_rho[gr](1,i0,j0,k0+1) +
                                   wghx*mMaterial_rho[gr](1,i0+1,j0,k0+1) ) +
   				   wghy*(    (1-wghx)*mMaterial_rho[gr](1,i0,j0+1,k0+1) + 
                                   wghx*mMaterial_rho[gr](1,i0+1,j0+1,k0+1) ) );

                       cp[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0) + 
                                   wghx*mMaterial_cp[gr](1,i0+1,j0,k0) ) +
                                   wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0) + 
                                   wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0) ) ) + 
                                   wghz*(  (1-wghy)*( (1-wghx)*mMaterial_cp[gr](1,i0,j0,k0+1) + 
                                   wghx*mMaterial_cp[gr](1,i0+1,j0,k0+1) ) +
                                   wghy*(     (1-wghx)*mMaterial_cp[gr](1,i0,j0+1,k0+1)+ 
                                   wghx*mMaterial_cp[gr](1,i0+1,j0+1,k0+1) ) );
       
                       cs[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0) + 
                                   wghx*mMaterial_cs[gr](1,i0+1,j0,k0) ) +
                                   wghy*(    (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0) + 
                                   wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0) ) ) + 
                                   wghz*(   (1-wghy)*( (1-wghx)*mMaterial_cs[gr](1,i0,j0,k0+1) + 
                                   wghx*mMaterial_cs[gr](1,i0+1,j0,k0+1) ) +
                                   wghy*(   (1-wghx)*mMaterial_cs[gr](1,i0,j0+1,k0+1)+ 
                                   wghx*mMaterial_cs[gr](1,i0+1,j0+1,k0+1) ) );

                       if( use_q ) {
                          xip[g](i, j, k) = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_qp[gr](1,i0+1,j0,k0) ) +
                                     wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*( (1-wghy)*(   (1-wghx)*mMaterial_qp[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_qp[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*( (1-wghx)*mMaterial_qp[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_qp[gr](1,i0+1,j0+1,k0+1) ) );

                          xis[g](i, j, k)= (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0) + 
                                     wghx*mMaterial_qs[gr](1,i0+1,j0,k0) ) +
                                     wghy*(      (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0) + 
                                     wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0) ) ) + 
                                     wghz*(   (1-wghy)*(   (1-wghx)*mMaterial_qs[gr](1,i0,j0,k0+1) + 
                                     wghx*mMaterial_qs[gr](1,i0+1,j0,k0+1) ) +
                                     wghy*(     (1-wghx)*mMaterial_qs[gr](1,i0,j0+1,k0+1)+ 
                                     wghx*mMaterial_qs[gr](1,i0+1,j0+1,k0+1) ) );
                       }

                       // Tang debug
                       /* float_sw4 rho_tmp = rho[g](i, j, k); */
                       /* ASSERT(rho[g](i, j, k) >= 1000); */
                       /* ASSERT(cp[g](i, j, k) >= 600); */
                       /* ASSERT(cs[g](i, j, k) >= 70); */
                       /* if (rho[g](i, j, k) > 3000 || rho[g](i, j, k) < 1500) */ 
                       /*     printf("Error with rho[%d]: %.2f\n", ind, rho[g](i, j, k)); */
                       /* if (cp[g](i, j, k) > 5700  || cp[g](i, j, k) < 600) */ 
                       /*     printf("Error with cp[%d]: %.2f\n", ind, cp[g](i, j, k)); */
                       /* if (cs[g](i, j, k) > 3400|| cs[g](i, j, k) < 70) */ 
                       /*     printf("Error with cs[%d]: %.2f\n", ind, cs[g](i, j, k)); */
                       /* if( use_q ) { */
                       /*   if (xip[g](i, j, k) > 377|| xip[g](i, j, k) < 12) */ 
                       /*       printf("Error with xip[%d]: %.2f\n", ind, xip[g](i, j, k)); */
                       /*   if (xis[ind] > 752|| xis[ind] < 25) */ 
                       /*       printf("Error with xis[%d]: %.2f\n", ind, xis[ind]); */
                       /* } */
     
                   } // End else
		} // End if inside
		else
		   outside++;
	    }

    /* for (int gk = mEW->m_kStart[g]; gk <= mEW->m_kEnd[g]; ++gk) */
    /* for (int gk = mEW->m_kStartInt[g]; gk <= mEW->m_kEndInt[g]; ++gk) */
    /*   for (int gj = mEW->m_jStartInt[g]; gj <= mEW->m_jEndInt[g]; ++gj) */
    /*     for (int gi = mEW->m_iStartInt[g]; gi <= mEW->m_iEndInt[g]; ++gi) { */
    /*       if (rho[g](gi,gj,gk) < 0) { */
    /*         printf("Error with grid %d rho(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, rho[g](gi,gj,gk), rho[g](gi,gj,gk+1)); */
    /*         rho[g](gi,gj,gk) = rho[g](gi,gj,gk+1); */
    /*       } */
    /*       if (cp[g](gi,gj,gk) < 100) { */
    /*         printf("Error with grid %d cp(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, cp[g](gi,gj,gk), cp[g](gi,gj,gk+1)); */
    /*         cp[g](gi,gj,gk) = cp[g](gi,gj,gk+1); */
    /*       } */
    /*       if (cs[g](gi,gj,gk) < 0) { */
    /*         printf("Error with grid %d cs(%d, %d, %d): %.2f, changing to %.2f\n", g, gi, gj, gk, cs[g](gi,gj,gk), cs[g](gi,gj,gk+1)); */
    /*         cs[g](gi,gj,gk) = cs[g](gi,gj,gk+1); */
    /*       } */
    /*       if( use_q ) { */
    /*         if (xip[g](gi,gj,gk) < 0) { */
    /*           printf("Error with grid %d xip(%d, %d, %d): %.2f, changing to %.2f\n", gi, g, gj, gk, xip[g](gi,gj,gk), xip[g](gi,gj,gk+1)); */
    /*           xip[g](gi,gj,gk) = xip[g](gi,gj,gk+1); */
    /*         } */
    /*         if (xis[g](gi,gj,gk) < 0) { */
    /*           printf("Error with grid %d xis(%d, %d, %d): %.2f, changing to %.2f\n", gi, g, gj, gk, xis[g](gi,gj,gk), xis[g](gi,gj,gk+1)); */
    /*           xis[g](gi,gj,gk) = xis[g](gi,gj,gk+1); */
    /*         } */
    /*       } */

    /*     } */

    if (mEW->getRank() == 0) 
      printf("After interpolation from sfile\n");

    printf("Rank %d grid %d: rho min = %.2f, max = %.2f\n", mEW->getRank(), g, rho[g].minimum(), rho[g].maximum());
    printf("Rank %d grid %d: cp min = %.2f, max = %.2f\n", mEW->getRank(), g, cp[g].minimum(), cp[g].maximum());
    printf("Rank %d grid %d: cs min = %.2f, max = %.2f\n", mEW->getRank(), g, cs[g].minimum(), cs[g].maximum());

    if( use_q ) {
      printf("Rank %d grid %d: xip min = %.2f, max = %.2f\n", mEW->getRank(), g, xip[g].minimum(), xip[g].maximum());
      printf("Rank %d grid %d: xis min = %.2f, max = %.2f\n", mEW->getRank(), g, xis[g].minimum(), xis[g].maximum());
    }

   } // end for g...

   mEW->communicate_arrays( rho );
   mEW->communicate_arrays( cs );
   mEW->communicate_arrays( cp );
   mEW->material_ic( rho );
   mEW->material_ic( cs );
   mEW->material_ic( cp );
   if( use_q ) {
      mEW->communicate_arrays( xis );
      mEW->communicate_arrays( xip );
      mEW->material_ic( xis);
      mEW->material_ic( xip );
   }

   size_t materialSum, outsideSum;
   int mpisizelong, mpisizelonglong, mpisizeint;
   MPI_Type_size(MPI_LONG,&mpisizelong );
   MPI_Type_size(MPI_LONG_LONG,&mpisizelonglong );
   MPI_Type_size(MPI_INT,&mpisizeint );
   if( sizeof(size_t) == mpisizelong ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizelonglong ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizeint ) {
      MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else {
      int materialsumi, outsidesumi, materiali=material, outsidei=outside;
      MPI_Reduce(&materiali, &materialsumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outsidei,   &outsidesumi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      materialSum=materialsumi;
      outsideSum=outsidesumi;
   }
   if (mEW->getRank() == 0)
      //      cout << endl 
      //           << "--------------------------------------------------------------\n"
      //           << "Sfile Initialized Node Types: " << endl
      //           << "   Material:        " << materialSum << endl
      //           << endl
      //           << "*Outside Domain:    " << outsideSum << endl
      //           << endl 
      //           << "--------------------------------------------------------------\n"
      //           << endl;
      cout << endl
	   << "sfile command: outside = " << outsideSum << ", material = " << materialSum << endl;
}


//-----------------------------------------------------------------------
void MaterialSfile::read_sfile2()
{
  // Timers
  double time_start = MPI_Wtime();

  string rname = "MaterialSfile::read_sfile2";

  // Figure out bounding box in this processor
  float_sw4 xmin=1e38, xmax=-1e38, ymin=1e38, ymax=-1e38, zmin=1e38, zmax=-1e38;
  for( int g=0 ; g < mEW->mNumberOfGrids ; g++ ) {
     float_sw4 h=mEW->mGridSize[g];
     if( xmin > (mEW->m_iStart[g]-1)*h )
        xmin =  (mEW->m_iStart[g]-1)*h;
     if( xmax < (mEW->m_iEnd[g]-1)*h )
        xmax =  (mEW->m_iEnd[g]-1)*h;
     if( ymin > (mEW->m_jStart[g]-1)*h )
        ymin =  (mEW->m_jStart[g]-1)*h;
     if( ymax < (mEW->m_jEnd[g]-1)*h )
        ymax =  (mEW->m_jEnd[g]-1)*h;
     if( mEW->topographyExists() && g == mEW->mNumberOfGrids-1 ) {
        int kb=mEW->m_kStart[g], ke=mEW->m_kEnd[g];
        for( int j=mEW->m_jStart[g] ; j <= mEW->m_jEnd[g] ; j++ ) {
           for( int i=mEW->m_iStart[g] ; i <= mEW->m_iEnd[g] ; i++ ) {
              if( zmin > mEW->mZ(i,j,kb) )
                  zmin = mEW->mZ(i,j,kb);
              if( zmax < mEW->mZ(i,j,ke) )
                  zmax = mEW->mZ(i,j,ke);
           }
        }
     }
     else {
        if( zmin > (mEW->m_kStart[g]-1)*h + mEW->m_zmin[g] ) 
           zmin = (mEW->m_kStart[g]-1)*h + mEW->m_zmin[g];
        if( zmax < (mEW->m_kEnd[g]-1)*h + mEW->m_zmin[g] )
           zmax = (mEW->m_kEnd[g]-1)*h + mEW->m_zmin[g];
     }
  }
  m_xminloc = xmin;
  m_xmaxloc = xmax;
  m_yminloc = ymin;
  m_ymaxloc = ymax;
  m_zminloc = zmin;
  m_zmaxloc = zmax;

  string fname = m_model_dir + "/" + m_model_file;

  hid_t file_id, dataset_id, datatype_id, h5_dtype, group_id, grid_id, memspace_id, filespace_id, attr_id, plist_id;
  int prec;
  double lonlataz[3];
  herr_t ierr;
  hsize_t dim[3];

  // Open file
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
  if (file_id < 0) {
     cout << "Could not open hdf5 file: " << fname.c_str()<< endl;
     MPI_Abort(MPI_COMM_WORLD, file_id);
  }

  // Read sfile header. Translate each patch into SW4 Cartesian coordinate system
  // Origin longitude, latitude, azimuth
  attr_id = H5Aopen(file_id, "Origin longitude, latitude, azimuth", H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, lonlataz);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

  // ---------- attenuation on file ?
  int att;
  attr_id = H5Aopen(file_id, "Attenuation", H5P_DEFAULT);
  ASSERT(attr_id >= 0);
  ierr = H5Aread(attr_id, H5T_NATIVE_INT, &att);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

  m_use_attenuation = (att==1);

  // ---------- azimuth on file
  double alpha = lonlataz[2], lon0 = lonlataz[0], lat0 = lonlataz[1];
  CHECK_INPUT( fabs(alpha-mEW->getGridAzimuth()) < 1e-6, "ERROR: sfile azimuth must be equal "
               "to coordinate system azimuth" << " azimuth on sfile = " << alpha << 
               " azimuth of coordinate sytem = " << mEW->getGridAzimuth() );

  // ---------- origin on file
  mEW->computeCartesianCoord( m_x0, m_y0, lon0, lat0 );

  // ---------- number of blocks on file
  attr_id = H5Aopen(file_id, "ngrids", H5P_DEFAULT);
  ierr = H5Aread(attr_id, H5T_NATIVE_INT, &m_npatches);
  ASSERT(ierr >= 0);
  H5Aclose(attr_id);

//test
  if (mEW->getRank()==0 && mEW->getVerbosity() >= 2) {
    printf("Sfile header: att=%i\n", att);
    printf("              azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
    printf("              nblocks=%i\n", m_npatches);
  }
     
  m_hh.resize(m_npatches);
  m_hv.resize(m_npatches);
  m_z0.resize(m_npatches);
  m_ni.resize(m_npatches);
  m_nj.resize(m_npatches);
  m_nk.resize(m_npatches);

  // ---------- read block headers
  vector<int> ncblock(m_npatches);
  char grid_name[16];

  group_id = H5Gopen(file_id, "Material_model", H5P_DEFAULT);
  ASSERT(group_id >= 0);
  for( int p=0 ; p < m_npatches ; p++ ) {

     sprintf(grid_name, "grid_%d", p);
     grid_id = H5Gopen(group_id, grid_name, H5P_DEFAULT);
     ASSERT(grid_id >= 0);

     // ---------- first part of block header
     double hs[3];
     attr_id = H5Aopen(grid_id, "Horizontal, vertical grid size, base z-level", H5P_DEFAULT);
     ASSERT(attr_id >= 0);
     ierr = H5Aread(attr_id, H5T_NATIVE_DOUBLE, hs);
     ASSERT(ierr >= 0);
     H5Aclose(attr_id);

     m_hh[p] = static_cast<float_sw4>(hs[0]);
     m_hv[p] = static_cast<float_sw4>(hs[1]);
     m_z0[p] = static_cast<float_sw4>(hs[2]);

     // ---------- second part of block header
     attr_id = H5Aopen(grid_id, "Number of components", H5P_DEFAULT);
     ASSERT(attr_id >= 0);
     ierr = H5Aread(attr_id, H5T_NATIVE_INT, &ncblock[p]);
     ASSERT(ierr >= 0);
     H5Aclose(attr_id);

     dataset_id = H5Dopen(grid_id, "Cp", H5P_DEFAULT);
     ASSERT(dataset_id >= 0);

     filespace_id = H5Dget_space(dataset_id);
     H5Sget_simple_extent_dims(filespace_id, dim, NULL);
     H5Sclose(filespace_id);
     H5Dclose(dataset_id);
     H5Gclose(grid_id);

     m_ni[p]    = (int)dim[0];
     m_nj[p]    = (int)dim[1];
     m_nk[p]    = (int)dim[2];
     if (mEW->getRank()==0 && mEW->getVerbosity() >= 2) {
       printf("  header block #%i\n", p);
       printf("  hh=%e, hv=%e, z0=%e\n", m_hh[p], m_hv[p], m_z0[p]);
       printf("  nc=%i, ni=%i, nj=%i, nk=%i\n", ncblock[p], m_ni[p], m_nj[p], m_nk[p]);
     }
  } // End for each patch

  // Intersect local grid with grid on sfile, assume all patches have the
  // same x- and y- extent. Assume patch=0 is topography, patches =1,..npatches-1
  // are ordered from top to bottom.
  float_sw4 xminrf = m_x0,    xmaxrf = m_x0+(m_ni[0]-1)*m_hh[0];
  float_sw4 yminrf = m_y0,    ymaxrf = m_y0+(m_nj[0]-1)*m_hh[0];
  float_sw4 zminrf = m_z0[0], zmaxrf = m_z0[m_npatches-1] + (m_nk[m_npatches-1]-1)*m_hv[m_npatches-1];

  if( xminrf > m_xminloc )
     m_xminloc = xminrf;
  if( xmaxrf < m_xmaxloc )
     m_xmaxloc = xmaxrf;
  if( yminrf > m_yminloc )
     m_yminloc = yminrf;
  if( ymaxrf < m_ymaxloc )
     m_ymaxloc = ymaxrf;
  if( zminrf > m_zminloc )
     m_zminloc = zminrf;
  if( zmaxrf < m_zmaxloc )
     m_zmaxloc = zmaxrf;
  
  
  /* mMaterial.resize(m_npatches); */
  mMaterial_rho.resize(m_npatches);
  mMaterial_cp.resize(m_npatches);
  mMaterial_cs.resize(m_npatches);
  mMaterial_qp.resize(m_npatches);
  mMaterial_qs.resize(m_npatches);

  m_ifirst.resize(m_npatches);
  m_ilast.resize(m_npatches);
  m_jfirst.resize(m_npatches);
  m_jlast.resize(m_npatches);
  m_kfirst.resize(m_npatches);
  m_klast.resize(m_npatches);

  /* size_t pos0 = 5*sizeof(int) + 3*sizeof(double)+ len*sizeof(char) + */
  /*    m_npatches*(3*sizeof(double)+4*sizeof(int)); */

  m_outside = m_xminloc >= m_xmaxloc || m_yminloc >= m_ymaxloc;
  m_isempty.resize(m_npatches);

  if( !m_outside ) {
     // each patch, figure out a box that encloses [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
     for( int p=0 ; p < m_npatches ; p++ ) {
        m_ifirst[p] = static_cast<int>(floor( 1 + (m_xminloc-m_x0)/m_hh[p]));
        m_ilast[p]  = static_cast<int>( ceil(  1 + (m_xmaxloc-m_x0)/m_hh[p]));
        m_jfirst[p] = static_cast<int>(floor( 1 + (m_yminloc-m_y0)/m_hh[p]));
        m_jlast[p]  = static_cast<int>( ceil(  1 + (m_ymaxloc-m_y0)/m_hh[p]));
        m_kfirst[p] = static_cast<int>(floor( 1 + (m_zminloc-m_z0[p])/m_hv[p]));
        m_klast[p]  = static_cast<int>( ceil(  1 + (m_zmaxloc-m_z0[p])/m_hv[p]));

        /* if( p == 0 ) */
           /* m_kfirst[0] = m_klast[0] = 1; //topography */

        // Limit index ranges to global size limits
        if( m_ifirst[p] < 1 )
           m_ifirst[p] = 1;
        if( m_ilast[p] > m_ni[p] )
           m_ilast[p] = m_ni[p];
        if( m_jfirst[p] < 1 )
           m_jfirst[p] = 1;
        if( m_jlast[p] > m_nj[p] )
           m_jlast[p] = m_nj[p];
        if( m_kfirst[p] < 1 )
           m_kfirst[p] = 1;
        if( m_klast[p] > m_nk[p] )
           m_klast[p] = m_nk[p];

        m_isempty[p] = false;
        if( m_klast[p] < m_kfirst[p] ) {
           m_isempty[p] = true;
           m_ilast[p] = 0;
           m_jlast[p] = 0;
           m_klast[p] = 0;
           m_ifirst[p] = 1;
           m_jfirst[p] = 1;
           m_kfirst[p] = 1;
        }
        if (mEW->getRank() == 0 && mEW->getVerbosity() >= 2) {
        /* if (mEW->getVerbosity() >= 3) { */
           cout << "myRank = " << mEW->getRank() << endl;
           cout << "patch nr " << p << " i " << m_ifirst[p] << " " << m_ilast[p] <<
    	  " j " << m_jfirst[p] << " " << m_jlast[p] << 
    	  " k " << m_kfirst[p] << " " << m_klast[p] << endl;
           cout << "nr components " << ncblock[p] << endl;
           cout << "patch nr global size " << m_ni[p] << " x " << m_nj[p] << " x " << m_nk[p] << endl;
        }
     }
  }
  else {
     for( int p=0 ; p < m_npatches ; p++ ) {
        m_isempty[p] = true;
        m_ilast[p] = 0;
        m_jlast[p] = 0;
        m_klast[p] = 0;
        m_ifirst[p] = 1;
        m_jfirst[p] = 1;
        m_kfirst[p] = 1;
     }
  }
  vector<int> isempty(m_npatches), isemptymin(m_npatches);
  for( int p=0 ; p < m_npatches ; p++ )
     isempty[p] = m_isempty[p];
  MPI_Allreduce( &isempty[0], &isemptymin[0], m_npatches, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  for( int p=0 ; p < m_npatches ; p++ )
     m_isempty[p] = (isemptymin[p] == 1);

  //Allocate memory
  for( int p=0 ; p < m_npatches ; p++ ) {
     try {
        if( !m_isempty[p] ) {
           /* mMaterial[p].define(ncblock[p],m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]); */
           mMaterial_rho[p].define(1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_cp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_cs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_qp[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
           mMaterial_qs[p].define( 1,m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
        }
     }
     catch( bad_alloc& ba ) {
        cout << "Processor " << mEW->getRank() << " allocation of mMaterial failed." << endl;
        cout << "p= "<< p << " ncblock= " << ncblock[p] << " ifirst,ilast " << m_ifirst[p] << " " << m_ilast[p] <<
           " jfirst,jlast " << m_jfirst[p] << " " << m_jlast[p] <<
           " kfirst,klast " << m_kfirst[p] << " " << m_klast[p] << 
           " Exception= " << ba.what() << endl;
        MPI_Abort(MPI_COMM_WORLD,0);
     }
  }

  // Read interfaces
  hid_t topo_grp, topo_id;
  mInterface.resize(m_npatches+1);
  char intf_name[32];
  void  *in_data;

  topo_grp = H5Gopen(file_id, "Z_interfaces", H5P_DEFAULT);
  ASSERT(topo_grp >= 0);

  for (int p = 0; p < m_npatches+1; p++) {
    sprintf(intf_name, "z_values_%d", p);
    dataset_id = H5Dopen(topo_grp, intf_name, H5P_DEFAULT);
    ASSERT(dataset_id >= 0);

    filespace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(filespace_id, dim, NULL);
    H5Sclose(filespace_id);

    mInterface[p].define(1, (int)dim[0], 1, (int)dim[1], 1, 1);

    float  *f_data = new  float[dim[0]*dim[1]];
    double *d_data = new double[dim[0]*dim[1]];

    if (p == 0) {
        // Get precision
        datatype_id = H5Dget_type(dataset_id);
        prec = (int)H5Tget_size(datatype_id);
        H5Tclose(datatype_id);
  
        if (prec == 4) 
            h5_dtype = H5T_NATIVE_FLOAT;
        else if (prec == 8) 
            h5_dtype = H5T_NATIVE_DOUBLE;
    }
    if (prec == 4) 
        in_data = (void*)f_data;
    else if (prec == 8) 
        in_data = (void*)d_data;
  
    ierr = H5Dread(dataset_id, h5_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, in_data);
    ASSERT(ierr >= 0);
    H5Dclose(dataset_id);

    if (prec == 4) { mInterface[p].assign(f_data); }
    else {           mInterface[p].assign(d_data); }

    printf("Rank %d: interface %d min = %.2f, max = %.2f\n", 
            mEW->getRank(), p,  mInterface[p].minimum(), mInterface[p].maximum());

    delete[] f_data;
    delete[] d_data;
  }
  H5Gclose(topo_grp);

  // Tang: TODO
  bool roworder = true;
  hid_t rho_id, cs_id, cp_id, qs_id, qp_id;
  for( int p = 0 ; p < m_npatches ; p++ ) {
    if( !m_isempty[p] ) {
      int global[3]={ m_ni[p], m_nj[p], m_nk[p] };
      int local[3] ={ m_ilast[p]-m_ifirst[p]+1, m_jlast[p]-m_jfirst[p]+1, m_klast[p]-m_kfirst[p]+1 };
      int start[3] ={ m_ifirst[p]-1, m_jfirst[p]-1, m_kfirst[p]-1 };
      /* if( roworder ) { */
      /*   int tmp=global[0]; */
      /*   global[0]=global[2]; */
      /*   global[2]=tmp; */
      /*   tmp=local[0]; */
      /*   local[0]=local[2]; */
      /*   local[2]=tmp; */
      /*   tmp=start[0]; */
      /*   start[0]=start[2]; */
      /*   start[2]=tmp; */
      /* } */

      printf("\nRank %d: start (%d, %d, %d), count (%d, %d, %d), global (%d, %d, %d), %d points\n\n",
              mEW->getRank(), start[0], start[1], start[2], local[0], local[1], local[2], 
              global[0], global[1], global[2], mMaterial_rho[p].m_npts);
      fflush(stdout);

      /* float_sw4 *material_dble = new float_sw4[mMaterial_rho[p].m_npts]; */
      float  *f_data = new float[mMaterial_rho[p].m_npts];
      double *d_data = new double[mMaterial_rho[p].m_npts];

      if (prec == 4) 
          in_data = (void*)f_data;
      else if (prec == 8) 
          in_data = (void*)d_data;


      sprintf(grid_name, "grid_%d", p);
      grid_id = H5Gopen(group_id, grid_name, H5P_DEFAULT);
      ASSERT(grid_id >= 0);

      rho_id = H5Dopen(grid_id, "Rho", H5P_DEFAULT);
      ASSERT(rho_id >= 0);
      cp_id  = H5Dopen(grid_id, "Cp", H5P_DEFAULT);
      ASSERT(cp_id >= 0);
      cs_id  = H5Dopen(grid_id, "Cs", H5P_DEFAULT);
      ASSERT(cs_id >= 0);
      if (m_use_attenuation) {
        qp_id  = H5Dopen(grid_id, "Qp", H5P_DEFAULT);
        ASSERT(qp_id >= 0);
        qs_id  = H5Dopen(grid_id, "Qs", H5P_DEFAULT);
        ASSERT(qs_id >= 0);
      }

      hsize_t h5_global[3], h5_count[3], h5_start[3];
      for (int i = 0; i < 3; i++) {
          h5_global[i] = (hsize_t)global[i];
          h5_count[i]  = (hsize_t)local[i];
          h5_start[i]  = (hsize_t)start[i];
      }

      memspace_id = H5Screate_simple(3, h5_count, NULL);
      ASSERT(memspace_id >= 0);

      filespace_id = H5Dget_space(rho_id);
      ASSERT(filespace_id >= 0);
      ierr = H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, h5_start, NULL, h5_count, NULL);
      ASSERT(ierr >= 0);

      printf("Rank %d patch %d, selected %lu points\n",  mEW->getRank(), p, H5Sget_select_npoints(filespace_id) );
 
      // Read each var individually 
      ierr = H5Dread(rho_id, h5_dtype, memspace_id, filespace_id, H5P_DEFAULT, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_rho[p].assign(f_data); }
      else {            mMaterial_rho[p].assign(d_data); }
        printf("Rank %d: rho min = %.2f, max = %.2f\n", 
                mEW->getRank(), mMaterial_rho[p].minimum(), mMaterial_rho[p].maximum());

      ierr = H5Dread(cp_id, h5_dtype, memspace_id, filespace_id, H5P_DEFAULT, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_cp[p].assign(f_data); }
      else {            mMaterial_cp[p].assign(d_data); }
        printf("Rank %d: cp min = %.2f, max = %.2f\n", 
                mEW->getRank(), mMaterial_cp[p].minimum(), mMaterial_cp[p].maximum());

      ierr = H5Dread(cs_id, h5_dtype, memspace_id, filespace_id, H5P_DEFAULT, in_data);
      ASSERT(ierr >= 0);
      if( prec == 4 ) { mMaterial_cs[p].assign(f_data); }
      else {            mMaterial_cs[p].assign(d_data); }
        printf("Rank %d: cs min = %.2f, max = %.2f\n", 
                mEW->getRank(), mMaterial_cs[p].minimum(), mMaterial_cs[p].maximum());

      if (m_use_attenuation) {
          ierr = H5Dread(qp_id, h5_dtype, memspace_id, filespace_id, H5P_DEFAULT, in_data);
          ASSERT(ierr >= 0);
          if( prec == 4 ) { mMaterial_qp[p].assign(f_data); }
          else {            mMaterial_qp[p].assign(d_data); }
          printf("Rank %d: qp min = %.2f, max = %.2f\n", 
                  mEW->getRank(), mMaterial_qp[p].minimum(), mMaterial_qp[p].maximum());

          ierr = H5Dread(qs_id, h5_dtype, memspace_id, filespace_id, H5P_DEFAULT, in_data);
          ASSERT(ierr >= 0);
          if( prec == 4 ) { mMaterial_qs[p].assign(f_data); }
          else {            mMaterial_qs[p].assign(d_data); }
          printf("Rank %d: qs min = %.2f, max = %.2f\n", 
                  mEW->getRank(), mMaterial_qs[p].minimum(), mMaterial_qs[p].maximum());
      }
      H5Sclose(memspace_id);
      H5Sclose(filespace_id);
      H5Dclose(rho_id);
      H5Dclose(cp_id);
      H5Dclose(cs_id);
      if (m_use_attenuation) {
          H5Dclose(qp_id);
          H5Dclose(qs_id);
      }


      /* mMaterial_rho[p].assign( material_dble, 0 ); */

      delete[] f_data;
      delete[] d_data;
      /* delete[] material_dble; */
      // Tang: is transpose needed?
      if( roworder ) {
         mMaterial_rho[p].transposeik();
         mMaterial_cp[p].transposeik();
         mMaterial_cs[p].transposeik();
         mMaterial_qp[p].transposeik();
         mMaterial_qs[p].transposeik();
      }

      H5Gclose(grid_id);
    } // End if !m_isempty
  } // End for

  H5Gclose(group_id);
  H5Fclose(file_id);

  fill_in_fluids();
  // material_check(false);

  double time_end = MPI_Wtime();
  if (mEW->getRank() == 0)
     cout << "MaterialSfile::read_sfile, time to read material file: " << time_end - time_start << " seconds." << endl;
  cout.flush();
}

/* //----------------------------------------------------------------------- */
/* void MaterialSfile::read_sfile() */
/* { */
/* #ifdef USE_HDF5 */
/*    string filename = m_model_dir + "/" + m_model_file; */
/*    SfileHDF5::read_sfile_material(filename, *mEW, *this, mMaterial, */ 
/*       mInterface); */

/* #else */
/* 	 cout << "WARNING: sw4 not compiled with hdf5=yes, " << */
/*     "--> ignoring MaterialSfile::read_sfile, no-op" << endl; */
/* #endif */
/* } */

/* //----------------------------------------------------------------------- */
/* void MaterialSfile::read_topo(const std::string &file, */ 
/*     const std::string &path, EW& ew, Sarray& gridElev, */
/*     float_sw4& lon0, float_sw4& lat0, float_sw4& azim, float_sw4& hh) */
/* { */
/* #ifdef USE_HDF5 */
/*    // Don't look in results path */
/*    // string filename = path + "/" + file; */
/*    SfileHDF5::read_sfile_topo(file, ew, gridElev, lon0, lat0, azim, hh); */
/* #else */
/* 	 cout << "WARNING: sw4 not compiled with hdf5=yes, " << */
/*     "--> ignoring MaterialSfile::read_topo, no-op" << endl; */
/* #endif */
/* } */

/* //----------------------------------------------------------------------- */
/* void MaterialSfile::write_sfile(const std::string &rfile_dir, */ 
/*     const std::string &rfile_name) */
/* { */
/* #ifdef USE_HDF5 */
/*    MaterialRfile rfile; */
/*    read_rfile(rfile, rfile_dir, rfile_name); */
/*    /1* */
/*    SfileHDF5::write_sfile(m_model_file, m_model_dir, */
/*       rfile, m_vec_depths, m_horizontalInterval); */
/*    *1/ */
/*    SfileHDF5::write_sfile(m_model_file, m_model_dir, */
/*       rfile.m_use_attenuation, rfile.mMaterial, */
/*       rfile.m_ni, rfile.m_nj, rfile.m_nk, */
/*       rfile.m_z0, rfile.m_hh, rfile.m_hv, */
/*       rfile.m_lon0, rfile.m_lat0, rfile.m_azim, */
/*       m_vec_depths, m_horizontalInterval); */
/* #else */
/* 	 cout << "WARNING: sw4 not compiled with hdf5=yes, " << */
/*     "--> ignoring MaterialSfile::write_sfile, no-op" << endl; */
/* #endif */
/* } */

/* //----------------------------------------------------------------------- */
/* void MaterialSfile::read_rfile(MaterialRfile& rfile, */ 
/*     const std::string &rfile_dir, const std::string &rfile_name) */
/* { */
/*   string rname = "MaterialSfile::read_rfile"; */

/*   // MaterialRfile is just a container for us to fill from rfile */
/*   rfile.m_model_dir = rfile_dir; */ 
/*   rfile.m_model_file = rfile_name; */ 
/*   rfile.read_whole_rfile(); */
/* } */

//-----------------------------------------------------------------------
void MaterialSfile::fill_in_fluids()
{
// Start from p=1, p=0 is the topography.
//   for( int p=1 ; p < m_npatches ; p++ )
// start from the last (bottom) block and progress upwards
   if( !m_outside )
   {
   for( int p=m_npatches-1 ; p >=1; p-- )
   {
      if( !m_isempty[p] )
      {
#pragma omp parallel for	 
	 for( int j=mMaterial_cs[p].m_jb ; j <= mMaterial_cs[p].m_je ; j++ )
	    for( int i=mMaterial_cs[p].m_ib ; i <= mMaterial_cs[p].m_ie ; i++ )
	    {
	       int k0 = mMaterial_cs[p].m_kb;
	       while( mMaterial_cs[p](1,i,j,k0) == -999 && k0 < mMaterial_cs[p].m_ke )
		  k0++;
// consider the case where the top block is all water. Then k0 = mMaterial[p].m_ke and mMaterial[p](3,i,j,k0)=-999
   // k0 is now the first k with cs > 0.
	       if (mMaterial_cs[p](1,i,j,k0)==-999)
	       {
// get value from block p+1
		  if (p<m_npatches-1)
		  {
		     int pd=p+1, id, jd, kd; // index of donor block
		     float_sw4 xm=(i-1)*m_hh[p];
		     float_sw4 ym=(j-1)*m_hh[p];
// get closest (id,jd) index on patch pd
		     id = static_cast<int>( 1 + trunc(xm/m_hh[pd]) );
		     jd = static_cast<int>( 1 + trunc(ym/m_hh[pd]) );
		     kd = mMaterial_cs[pd].m_kb; // get value from top of block pd
		
		     if (! (id >= mMaterial_cs[pd].m_ib && id <= mMaterial_cs[pd].m_ie && 
			    jd >= mMaterial_cs[pd].m_jb && jd <= mMaterial_cs[pd].m_je ))
		     {
// out of bounds: find nearest interior point
			if (id < mMaterial_cs[pd].m_ib) id=mMaterial_cs[pd].m_ib;
			if (id > mMaterial_cs[pd].m_ie) id=mMaterial_cs[pd].m_ie;
			if (jd < mMaterial_cs[pd].m_jb) jd=mMaterial_cs[pd].m_jb;
			if (jd > mMaterial_cs[pd].m_je) jd=mMaterial_cs[pd].m_je;
			
			printf("WARNING: nearest grid point to (%e,%e) was outside local part of block pd=%i\n"
			 " using id=%i, jd=%i, at (%e, %e)\n", xm, ym, pd, id, jd, (id-1)*m_hh[pd], (jd-1)*m_hh[pd]);

		     }
// get values from block 'pd'
		     mMaterial_rho[p](1,i,j,k0)= mMaterial_rho[pd](1,id,jd,kd);
		     mMaterial_cp[p](1,i,j,k0)= mMaterial_cp[pd](1,id,jd,kd);
		     mMaterial_cs[p](1,i,j,k0)= mMaterial_cs[pd](1,id,jd,kd);
		     if (m_use_attenuation)
		     {
			mMaterial_qp[p](1,i,j,k0)= mMaterial_qp[pd](1,id,jd,kd);
			mMaterial_qs[p](1,i,j,k0)= mMaterial_qs[pd](1,id,jd,kd);
		     }
		  }
		  else
		  {
		     printf("ERROR: found undefined material properties in last material block\n"
			    " patch p=%i, i=%i, j=%i, k0=%i\n", p, i, j, k0);
		  }
	       }
	       for( int k=mMaterial_cs[p].m_kb ; k < k0 ; k++ )
	       {
		  mMaterial_rho[p](1,i,j,k) = mMaterial_rho[p](1,i,j,k0);
		  mMaterial_cp[p](1,i,j,k)  = mMaterial_cp[p](1,i,j,k0);
		  mMaterial_cs[p](1,i,j,k)  = mMaterial_cs[p](1,i,j,k0);
		  if( m_use_attenuation )
		  {
		     mMaterial_qp[p](1,i,j,k) = mMaterial_qp[p](1,i,j,k0);
		     mMaterial_qs[p](1,i,j,k) = mMaterial_qs[p](1,i,j,k0);
		  }
	       }
	    }
      }
   }
   }   
}

//-----------------------------------------------------------------------
void MaterialSfile::material_check( bool water )
{
   bool printsmallcpcs=false;
   for( int p=1 ; p < m_npatches ; p++ )
   {
      double csmin=1e38,cpmin=1e38,cratmin=1e38,csmax=-1e38,cpmax=-1e38,cratmax=-1e38;
      double rhomin=1e38, rhomax=-1e38;
      for( int k=mMaterial_cs[p].m_kb ; k<= mMaterial_cs[p].m_ke ; k++ )
	 for( int j=mMaterial_cs[p].m_jb ; j<= mMaterial_cs[p].m_je ; j++ )
	    for( int i=mMaterial_cs[p].m_ib ; i<= mMaterial_cs[p].m_ie ; i++ )
	    {
	       if( water || mMaterial_cs[p](1,i,j,k) != -999 )
	       {
		  if( mMaterial_rho[p](1,i,j,k) < rhomin )
		     rhomin = mMaterial_rho[p](1,i,j,k);
		  if( mMaterial_rho[p](1,i,j,k) > rhomax )
		     rhomax = mMaterial_rho[p](1,i,j,k);
		  if( mMaterial_cs[p](1,i,j,k) < csmin )
		     csmin = mMaterial_cs[p](1,i,j,k);
		  if( mMaterial_cs[p](1,i,j,k) > csmax )
		     csmax = mMaterial_cs[p](1,i,j,k);
		  if( mMaterial_cp[p](1,i,j,k) < cpmin )
		     cpmin = mMaterial_cp[p](1,i,j,k);
		  if( mMaterial_cp[p](1,i,j,k) > cpmax )
		     cpmax = mMaterial_cp[p](1,i,j,k);
		  double crat = mMaterial_cp[p](1,i,j,k)/mMaterial_cs[p](1,i,j,k);
		  if( crat < cratmin ) {
		     cratmin = crat;
		     if( printsmallcpcs && crat < 1.41 ) {
			cout << "crat= " << crat << " at " << i << " " <<  j << " " << k << endl;
			cout << " material is " << mMaterial_rho[p](1,i,j,k) << " " << mMaterial_cp[p](1,i,j,k) << " " 
			     << mMaterial_cs[p](1,i,j,k) << " " << mMaterial_qp[p](1,i,j,k) << " " << mMaterial_qs[p](1,i,j,k) << endl;
		     }
		  }
		  if( crat > cratmax )
		     cratmax = crat;
	       }
	    }
      double cmins[4]={csmin,cpmin,cratmin,rhomin}, cmaxs[4]={csmax,cpmax,cratmax,rhomax};
      double cminstot[4], cmaxstot[4];
      MPI_Reduce(cmins, cminstot, 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
      MPI_Reduce(cmaxs, cmaxstot, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
      int myid;
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      if( myid == 0 )
	 //	 if( mEW->getRank()==0 )
      {
	 if( p== 1 && !water )
	    cout << "R-file limits, away from water: " << endl;
	 else if( p== 1 )
	    cout << "R-file limits : " << endl;
	 cout << "  Patch no " << p << " : " << endl;
	 cout << "    cp    min and max " << cminstot[1] << " " << cmaxstot[1] << endl;
	 cout << "    cs    min and max " << cminstot[0] << " " << cmaxstot[0] << endl;
	 cout << "    cp/cs min and max " << cminstot[2] << " " << cmaxstot[2] << endl;
	 cout << "    rho   min and max " << cminstot[3] << " " << cmaxstot[3] << endl;
      }
   }
}


