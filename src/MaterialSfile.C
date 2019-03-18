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
void MaterialSfile::set_material_properties(std::vector<Sarray> & rho, 
			       std::vector<Sarray> & cs,
			       std::vector<Sarray> & cp, 
			       std::vector<Sarray> & xis, 
			       std::vector<Sarray> & xip )
{

// Assume attenuation arrays defined on all grids if they are defined on grid zero.
   bool use_q = m_use_attenuation && xis[0].is_defined() && xip[0].is_defined();
   size_t outside=0, material=0;
   for( int g=0 ; g < mEW->mNumberOfGrids  ; g++ )
   {
      bool curvilinear = mEW->topographyExists() && g == mEW->mNumberOfGrids-1;
      float_sw4* rhop=rho[g].c_ptr();
      float_sw4*  csp=cs[g].c_ptr();
      float_sw4*  cpp=cp[g].c_ptr();
      float_sw4* qp, *qs;
      if( use_q  )
      {
         qs = xis[g].c_ptr();
         qp = xip[g].c_ptr();
      }
      size_t ni=mEW->m_iEnd[g]-mEW->m_iStart[g]+1;
      size_t nj=mEW->m_jEnd[g]-mEW->m_jStart[g]+1;
      ssize_t ofs = -mEW->m_iStart[g]-ni*(mEW->m_jStart[g])-ni*nj*(mEW->m_kStart[g]);
#pragma omp parallel for reduction(+:material,outside)
      for (int k = mEW->m_kStart[g]; k <= mEW->m_kEnd[g]; ++k)
	 for (int j = mEW->m_jStartInt[g]; j <= mEW->m_jEndInt[g]; ++j)
	    for (int i = mEW->m_iStartInt[g]; i <= mEW->m_iEndInt[g]; ++i)
	    {
		float_sw4 x = (i-1)*mEW->mGridSize[g];
		float_sw4 y = (j-1)*mEW->mGridSize[g];
		float_sw4 z;
		if( curvilinear )
		   z = mEW->mZ(i,j,k);
		else
		   z = mEW->m_zmin[g] + (k-1)*mEW->mGridSize[g];
		size_t ind = ofs + i + ni*j + ni*nj*k;
		if( inside( x, y, z )  )
		{
		   material++;
                   int gr = m_npatches-1;
		   while( gr >= 1 && z < m_z0[gr] )
		      gr--;
                   if( gr==0 ) // gr=0 is topography.
		      gr = 1;
		   int i0 = static_cast<int>( trunc( 1 + (x-m_x0)/m_hh[gr]     ) );
		   int j0 = static_cast<int>( trunc( 1 + (y-m_y0)/m_hh[gr]     ) );
		   int k0 = static_cast<int>( trunc( 1 + (z-m_z0[gr])/m_hv[gr] ) ); 

		   // Use bilinear interpolation always:
                   bool intp_cubic = false;

	   // Bias stencil near the boundary, need to communicate arrays afterwards.
		   if( i0 <= m_ifirst[gr] )
		   {
		      i0 = m_ifirst[gr];
                      intp_cubic = false;
		   }
		   if( i0 >= m_ilast[gr]-1 )
		   {
		      i0 = m_ilast[gr]-1;
                      intp_cubic = false;
		   }
		   if( j0 <= m_jfirst[gr] )
		   {
		      j0 = m_jfirst[gr];
		      intp_cubic = false;
		   }
		   if( j0 >= m_jlast[gr]-1 )
		   {
		      j0 = m_jlast[gr]-1;
                      intp_cubic = false;
		   }
		   if( k0 <= m_kfirst[gr] )
		   {
		      k0 = m_kfirst[gr];
                      intp_cubic = false;
		   }
		   if( k0 >= m_klast[gr]-1 )
		   {
		      k0 = m_klast[gr]-1;
                      intp_cubic = false;
		   }

       ASSERT(!intp_cubic);
       // bilinear intp.
       float_sw4 wghx = (x-( (i0-1)*m_hh[gr]+m_x0) )/m_hh[gr];
       float_sw4 wghy = (y-( (j0-1)*m_hh[gr]+m_y0) )/m_hh[gr];
       float_sw4 wghz = (z-( (k0-1)*m_hv[gr]+m_z0[gr]) )/m_hv[gr];
       rhop[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial[gr](1,i0,j0,k0)    +wghx*mMaterial[gr](1,i0+1,j0,k0)   ) +
            wghy*(     (1-wghx)*mMaterial[gr](1,i0,j0+1,k0)  +wghx*mMaterial[gr](1,i0+1,j0+1,k0) ) ) + 
                       wghz*( (1-wghy)*( (1-wghx)*mMaterial[gr](1,i0,j0,k0+1)  +wghx*mMaterial[gr](1,i0+1,j0,k0+1) ) +
             wghy*(    (1-wghx)*mMaterial[gr](1,i0,j0+1,k0+1)+wghx*mMaterial[gr](1,i0+1,j0+1,k0+1) ) );

       cpp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial[gr](2,i0,j0,k0)    + wghx*mMaterial[gr](2,i0+1,j0,k0) ) +
           wghy*(     (1-wghx)*mMaterial[gr](2,i0,j0+1,k0)  + wghx*mMaterial[gr](2,i0+1,j0+1,k0) ) ) + 
                     wghz*(  (1-wghy)*( (1-wghx)*mMaterial[gr](2,i0,j0,k0+1)  + wghx*mMaterial[gr](2,i0+1,j0,k0+1) ) +
           wghy*(     (1-wghx)*mMaterial[gr](2,i0,j0+1,k0+1)+ wghx*mMaterial[gr](2,i0+1,j0+1,k0+1) ) );
 
       csp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial[gr](3,i0,j0,k0)    + wghx*mMaterial[gr](3,i0+1,j0,k0) ) +
            wghy*(    (1-wghx)*mMaterial[gr](3,i0,j0+1,k0)  + wghx*mMaterial[gr](3,i0+1,j0+1,k0) ) ) + 
                    wghz*(   (1-wghy)*( (1-wghx)*mMaterial[gr](3,i0,j0,k0+1)  + wghx*mMaterial[gr](3,i0+1,j0,k0+1) ) +
             wghy*(   (1-wghx)*mMaterial[gr](3,i0,j0+1,k0+1)+ wghx*mMaterial[gr](3,i0+1,j0+1,k0+1) ) );
       if( use_q )
       {
          qp[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial[gr](4,i0,j0,k0)    + wghx*mMaterial[gr](4,i0+1,j0,k0) ) +
            wghy*(      (1-wghx)*mMaterial[gr](4,i0,j0+1,k0)  + wghx*mMaterial[gr](4,i0+1,j0+1,k0) ) ) + 
                    wghz*(   (1-wghy)*(   (1-wghx)*mMaterial[gr](4,i0,j0,k0+1)  + wghx*mMaterial[gr](4,i0+1,j0,k0+1) ) +
             wghy*(     (1-wghx)*mMaterial[gr](4,i0,j0+1,k0+1)+ wghx*mMaterial[gr](4,i0+1,j0+1,k0+1) ) );
          qs[ind] = (1-wghz)*( (1-wghy)*( (1-wghx)*mMaterial[gr](5,i0,j0,k0)    + wghx*mMaterial[gr](5,i0+1,j0,k0) ) +
            wghy*(      (1-wghx)*mMaterial[gr](5,i0,j0+1,k0)  + wghx*mMaterial[gr](5,i0+1,j0+1,k0) ) ) + 
                    wghz*(   (1-wghy)*(   (1-wghx)*mMaterial[gr](5,i0,j0,k0+1)  + wghx*mMaterial[gr](5,i0+1,j0,k0+1) ) +
             wghy*(     (1-wghx)*mMaterial[gr](5,i0,j0+1,k0+1)+ wghx*mMaterial[gr](5,i0+1,j0+1,k0+1) ) );
       }
       } // if inside
     } // loop over points
   } // end for g...
#if 0
   mEW->communicate_arrays( rho );
   mEW->communicate_arrays( cs );
   mEW->communicate_arrays( cp );
   mEW->material_ic( rho );
   mEW->material_ic( cs );
   mEW->material_ic( cp );
   if( use_q )
   {
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
   if( sizeof(size_t) == mpisizelong )
   {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizelonglong )
   {
      MPI_Reduce(&material, &materialSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else if( sizeof(size_t) == mpisizeint )
   {
      MPI_Reduce(&material, &materialSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce(&outside,   &outsideSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
   }
   else
   {
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
      
#endif // #if 0
}


//-----------------------------------------------------------------------
void MaterialSfile::read_sfile()
{
   string rname = "MaterialSfile::read_sfile";

   // Read everything from the sfile
   read_materials();

#if 0

      // ---------- length of projection string
      int len;
      nr = read( fd, &len, sizeof(int));
      if( nr != sizeof(int) )
      {
	 cout << rname << " Error reading len, nr= " << nr
	      << "bytes read" << endl;
         close(fd);
	 return;
      }
      if( swapbytes )
	 bswap.byte_rev( &len, 1, "int" );

      // ---------- skip projection string
      nr = lseek(fd, len*sizeof(char), SEEK_CUR );

      // ---------- number of blocks on file
      nr = read( fd, &m_npatches, sizeof(int) );
      if( nr != sizeof(int) )
      {
	 cout << rname << " Error reading npatches, nr= " << nr
	      << "bytes read" << endl;
         close(fd);
	 return;
      }
      if( swapbytes )
	 bswap.byte_rev( &m_npatches, 1, "int" );

// test
      if (mEW->getRank()==0 && mEW->getVerbosity() >= 2)
      {
	printf("Sfile header: magic=%i, prec=%i, att=%i\n", magic, prec, att);
	printf("              azimuth=%e, lon0=%e, lat0=%e\n", alpha, lon0, lat0);
	printf("              pstring-len=%i, pstr='%s'\n", len, "not implemented");
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
      for( int p=0 ; p < m_npatches ; p++ )
      {
	 // ---------- first part of block header
	 double hs[3];
	 nr = read( fd, hs, 3*sizeof(double) );
	 if( nr != 3*sizeof(double) )
	 {
	    cout << rname << " Error reading block spacings nr= " << nr
		 << "bytes read" << endl;
	    close(fd);
	    return;
	 }
	 if( swapbytes )
	    bswap.byte_rev( hs, 3, "double" );
	 m_hh[p] = static_cast<float_sw4>(hs[0]);
	 m_hv[p] = static_cast<float_sw4>(hs[1]);
	 m_z0[p] = static_cast<float_sw4>(hs[2]);

      // ---------- second part of block header
	 int dim[4];
	 nr = read( fd, dim, 4*sizeof(int) );
	 if( nr != 4*sizeof(int) )
	 {
	    cout << rname << " Error reading block dimensions nr= " << nr
		 << "bytes read" << endl;
	    close(fd);
	    return;
	 }
	 if( swapbytes )
	    bswap.byte_rev( dim, 4, "int" );

	 ncblock[p] = dim[0];
	 m_ni[p]    = dim[1];
	 m_nj[p]    = dim[2];
	 m_nk[p]    = dim[3];
// test
	 if (mEW->getRank()==0 && mEW->getVerbosity() >= 2)
	 {
	   printf("  header block #%i\n", p);
	   printf("  hh=%e, hv=%e, z0=%e\n", m_hh[p], m_hv[p], m_z0[p]);
	   printf("  nc=%i, ni=%i, nj=%i, nk=%i\n", ncblock[p], m_ni[p], m_nj[p], m_nk[p]);
	 }
      }

      // Intersect local grid with grid on sfile, assume all patches have the
      // same x- and y- extent. Assume patch=0 is topography, patches =1,..npatches-1
      // are ordered from top to bottom.
      float_sw4 xminrf = m_x0,    xmaxrf = m_x0+(m_ni[0]-1)*m_hh[0];
      float_sw4 yminrf = m_y0,    ymaxrf = m_y0+(m_nj[0]-1)*m_hh[0];
      float_sw4 zminrf = m_z0[1], zmaxrf = m_z0[m_npatches-1] + (m_nk[m_npatches-1]-1)*m_hv[m_npatches-1];

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
      
      
      mMaterial.resize(m_npatches);
      m_ifirst.resize(m_npatches);
      m_ilast.resize(m_npatches);
      m_jfirst.resize(m_npatches);
      m_jlast.resize(m_npatches);
      m_kfirst.resize(m_npatches);
      m_klast.resize(m_npatches);

      size_t pos0 = 5*sizeof(int) + 3*sizeof(double)+ len*sizeof(char) +
	 m_npatches*(3*sizeof(double)+4*sizeof(int));

      m_outside = m_xminloc >= m_xmaxloc || m_yminloc >= m_ymaxloc;
      m_isempty.resize(m_npatches);

      if( !m_outside )
      {
// For each patch, figure out a box that encloses [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
	 for( int p=0 ; p < m_npatches ; p++ )
	 {
	    m_ifirst[p] = static_cast<int>(floor( 1 + (m_xminloc-m_x0)/m_hh[p]));
	    m_ilast[p]  = static_cast<int>( ceil(  1 + (m_xmaxloc-m_x0)/m_hh[p]));
	    m_jfirst[p] = static_cast<int>(floor( 1 + (m_yminloc-m_y0)/m_hh[p]));
	    m_jlast[p]  = static_cast<int>( ceil(  1 + (m_ymaxloc-m_y0)/m_hh[p]));
	    m_kfirst[p] = static_cast<int>(floor( 1 + (m_zminloc-m_z0[p])/m_hv[p]));
	    m_klast[p]  = static_cast<int>( ceil(  1 + (m_zmaxloc-m_z0[p])/m_hv[p]));

	    if( p == 0 )
	       m_kfirst[0] = m_klast[0] = 1; //topography

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
	    if( m_klast[p] < m_kfirst[p] )
	    {
	       m_isempty[p] = true;
	       m_ilast[p] = 0;
	       m_jlast[p] = 0;
	       m_klast[p] = 0;
	       m_ifirst[p] = 1;
	       m_jfirst[p] = 1;
	       m_kfirst[p] = 1;
	    }
	    if (mEW->getVerbosity() >=3)
	    {
	       cout << "myRank = " << mEW->getRank() << endl;
	       cout << "patch nr " << p << " i " << m_ifirst[p] << " " << m_ilast[p] <<
		  " j " << m_jfirst[p] << " " << m_jlast[p] << 
		  " k " << m_kfirst[p] << " " << m_klast[p] << endl;
	       cout << "nr components " << ncblock[p] << endl;
	       cout << "patch nr global size " << m_ni[p] << " x " << m_nj[p] << " x " << m_nk[p] << endl;
	    }
	 
	 }
      }
      else
      {
	 for( int p=0 ; p < m_npatches ; p++ )
	 {
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

// Allocate memory
      for( int p=0 ; p < m_npatches ; p++ )
      {
	 try {
	    if( !m_isempty[p] )
	       mMaterial[p].define(ncblock[p],m_ifirst[p],m_ilast[p],m_jfirst[p],m_jlast[p],m_kfirst[p],m_klast[p]);
	 }
	 catch( bad_alloc& ba )
	 {
	    cout << "Processor " << mEW->getRank() << " allocation of mMaterial failed." << endl;
	    cout << "p= "<< p << " ncblock= " << ncblock[p] << " ifirst,ilast " << m_ifirst[p] << " " << m_ilast[p] <<
	       " jfirst,jlast " << m_jfirst[p] << " " << m_jlast[p] <<
	       " kfirst,klast " << m_kfirst[p] << " " << m_klast[p] << 
	       " Exception= " << ba.what() << endl;
	    MPI_Abort(MPI_COMM_WORLD,0);
	    
	 }
      }
      int iread = io_processor();

      //      if( mEW->getRank() == 7423 )
      //      {
      //	 cout << "DEBUG: from 7423, iread= " << iread << endl;
      //	 for( int p=0 ; p < m_npatches ; p++ )
      //	 {
      //	    cout << "p= "<< p << " ncblock= " << ncblock[p] << " ifirst,ilast " << m_ifirst[p] << " " << m_ilast[p] <<
      //	       " jfirst,jlast " << m_jfirst[p] << " " << m_jlast[p] <<
      //	       " kfirst,klast " << m_kfirst[p] << " " << m_klast[p] << endl;
      //	 }
      //      }

      //      vector<Parallel_IO*> pio(m_npatches);
      //      int bufsize =  200000;
      bool roworder = true;
      for( int p=0 ; p < m_npatches ; p++ )
      {
	 if( !m_isempty[p] )
	 {
	    int global[3]={ m_ni[p], m_nj[p], m_nk[p] };
	    int local[3] ={ m_ilast[p]-m_ifirst[p]+1, m_jlast[p]-m_jfirst[p]+1, m_klast[p]-m_kfirst[p]+1 };
	    int start[3] ={ m_ifirst[p]-1, m_jfirst[p]-1, m_kfirst[p]-1 };
            if( roworder )
	    {
	       int tmp=global[0];
	       global[0]=global[2];
	       global[2]=tmp;
	       tmp=local[0];
	       local[0]=local[2];
	       local[2]=tmp;
	       tmp=start[0];
	       start[0]=start[2];
	       start[2]=tmp;
	    }
	    Parallel_IO* pio = new Parallel_IO( iread, mEW->usingParallelFS(), global, local, start, m_bufsize );
	 //	 pio[p] = new Parallel_IO( iread, mEW->usingParallelFS(), global, local, start );
// Read corresponding part of patches
	    float_sw4* material_dble = new float_sw4[mMaterial[p].m_npts];
	    if( prec == 8 )
	       pio->read_array( &fd, ncblock[p], material_dble, pos0, "double", swapbytes );
//	       pio->read_array( &fd, ncblock[p], mMaterial[p].c_ptr(), pos0, "double", swapbytes );
	    else
	       pio->read_array( &fd, ncblock[p], material_dble, pos0, "float", swapbytes );
//	       pio->read_array( &fd, ncblock[p], mMaterial[p].c_ptr(), pos0, "float", swapbytes );
	    delete pio;
	    mMaterial[p].assign( material_dble, 0 );
	    delete[] material_dble;
	    if( roworder )
	       mMaterial[p].transposeik();
	 }
	 pos0 += ncblock[p]*m_ni[p]*static_cast<size_t>(m_nj[p])*m_nk[p]*flsize;
      }
      close(fd);

      fill_in_fluids();
   //      material_check(false);

      //      for( int p=0 ; p < m_npatches ; p++ )
      //	 delete pio[p];
      
      //      cout << "origin " << m_x0 << " " << m_y0 << endl;
      //      cout << "DEBUG -----" << endl;
      //      cout << "Block nr 2  " << endl;
      //      cout << "index " << 1+(50000-m_x0)/m_hh[1] << endl;
      //
      //      int i=26;
      //      for( int k=m_kfirst[1] ; k <= m_klast[1] ; k++ )
      //	 for( int j=m_jfirst[1] ; j <= m_jlast[1] ; j++ )
      //	    cout << j << " " << k << " " << mMaterial[1](1,i,j,k) << endl;
      //
      //      cout << "END DEBUG -----" << endl;
   }
   else
      cout << "MaterialSfile::read_sfile, error could not open file " << fname << endl;
#endif
}

#if 0
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
	 for( int j=mMaterial[p].m_jb ; j <= mMaterial[p].m_je ; j++ )
	    for( int i=mMaterial[p].m_ib ; i <= mMaterial[p].m_ie ; i++ )
	    {
	       int k0 = mMaterial[p].m_kb;
	       while( mMaterial[p](3,i,j,k0) == -999 && k0 < mMaterial[p].m_ke )
		  k0++;
// consider the case where the top block is all water. Then k0 = mMaterial[p].m_ke and mMaterial[p](3,i,j,k0)=-999
   // k0 is now the first k with cs > 0.
	       if (mMaterial[p](3,i,j,k0)==-999)
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
		     kd = mMaterial[pd].m_kb; // get value from top of block pd
		
		     if (! (id >= mMaterial[pd].m_ib && id <= mMaterial[pd].m_ie && 
			    jd >= mMaterial[pd].m_jb && jd <= mMaterial[pd].m_je ))
		     {
// out of bounds: find nearest interior point
			if (id < mMaterial[pd].m_ib) id=mMaterial[pd].m_ib;
			if (id > mMaterial[pd].m_ie) id=mMaterial[pd].m_ie;
			if (jd < mMaterial[pd].m_jb) jd=mMaterial[pd].m_jb;
			if (jd > mMaterial[pd].m_je) jd=mMaterial[pd].m_je;
			
			printf("WARNING: nearest grid point to (%e,%e) was outside local part of block pd=%i\n"
			 " using id=%i, jd=%i, at (%e, %e)\n", xm, ym, pd, id, jd, (id-1)*m_hh[pd], (jd-1)*m_hh[pd]);

		     }
// get values from block 'pd'
		     mMaterial[p](1,i,j,k0)= mMaterial[pd](1,id,jd,kd);
		     mMaterial[p](2,i,j,k0)= mMaterial[pd](2,id,jd,kd);
		     mMaterial[p](3,i,j,k0)= mMaterial[pd](3,id,jd,kd);
		     if (m_use_attenuation)
		     {
			mMaterial[p](4,i,j,k0)= mMaterial[pd](4,id,jd,kd);
			mMaterial[p](5,i,j,k0)= mMaterial[pd](5,id,jd,kd);
		     }
		  }
		  else
		  {
		     printf("ERROR: found undefined material properties in last material block\n"
			    " patch p=%i, i=%i, j=%i, k0=%i\n", p, i, j, k0);
		  }
	       }
	       for( int k=mMaterial[p].m_kb ; k < k0 ; k++ )
	       {
		  mMaterial[p](1,i,j,k) = mMaterial[p](1,i,j,k0);
		  mMaterial[p](2,i,j,k) = mMaterial[p](2,i,j,k0);
		  mMaterial[p](3,i,j,k) = mMaterial[p](3,i,j,k0);
		  if( mMaterial[p].m_nc > 3 )
		  {
		     mMaterial[p](4,i,j,k) = mMaterial[p](4,i,j,k0);
		     mMaterial[p](5,i,j,k) = mMaterial[p](5,i,j,k0);
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
      for( int k=mMaterial[p].m_kb ; k<= mMaterial[p].m_ke ; k++ )
	 for( int j=mMaterial[p].m_jb ; j<= mMaterial[p].m_je ; j++ )
	    for( int i=mMaterial[p].m_ib ; i<= mMaterial[p].m_ie ; i++ )
	    {
	       if( water || mMaterial[p](3,i,j,k) != -999 )
	       {
		  if( mMaterial[p](1,i,j,k) < rhomin )
		     rhomin = mMaterial[p](1,i,j,k);
		  if( mMaterial[p](1,i,j,k) > rhomax )
		     rhomax = mMaterial[p](1,i,j,k);
		  if( mMaterial[p](3,i,j,k) < csmin )
		     csmin = mMaterial[p](3,i,j,k);
		  if( mMaterial[p](3,i,j,k) > csmax )
		     csmax = mMaterial[p](3,i,j,k);
		  if( mMaterial[p](2,i,j,k) < cpmin )
		     cpmin = mMaterial[p](2,i,j,k);
		  if( mMaterial[p](2,i,j,k) > cpmax )
		     cpmax = mMaterial[p](2,i,j,k);
		  double crat = mMaterial[p](2,i,j,k)/mMaterial[p](3,i,j,k);
		  if( crat < cratmin )
		  {
		     cratmin = crat;
		     if( printsmallcpcs && crat < 1.41 )
		     {
			cout << "crat= " << crat << " at " << i << " " <<  j << " " << k << endl;
			cout << " material is " << mMaterial[p](1,i,j,k) << " " << mMaterial[p](2,i,j,k) << " " 
			     << mMaterial[p](3,i,j,k) << " " << mMaterial[p](4,i,j,k) << " " << mMaterial[p](5,i,j,k) << endl;
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
#endif

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
void MaterialSfile::read_materials()
{
#ifdef USE_HDF5
   string filename = m_model_dir + "/" + m_model_file;
   SfileHDF5::read_sfile_material(filename, *mEW, *this, mMaterial, 
      mInterface);

#else
	 cout << "WARNING: sw4 not compiled with hdf5=yes, " <<
    "--> ignoring MaterialSfile::read_materials, no-op" << endl;
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
