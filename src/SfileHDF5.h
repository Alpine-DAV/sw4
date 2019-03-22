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
#ifndef SW4_SFILEHDF5_H
#define SW4_SFILEHDF5_H

#include <string>
#include <sstream>
#include <vector>

#include "EW.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

class SfileHDF5
{
public:
  static void write_sfile(const std::string& file, const std::string& path, 
      EW& ew, MaterialRfile& model, vector<float_sw4>& depths, 
      int horizontalInterval=1);
  static void read_sfile_topo(const std::string& file, EW& ew, Sarray& gridElev,
      float_sw4& lon0, float_sw4& lat0, float_sw4& azim, float_sw4& hh);
  static void read_sfile_material(const std::string& file, EW& ew, 
      MaterialSfile& model,  vector<Sarray>& materials,
      vector<Sarray>& interfaces);

  // Structure for header data broadcast
  struct sfile_hdr
  {
    float h;
    float lon0;
    float lat0;
    float azim;
    int ngrids;
  };

protected:
  // Structure for keeping track of depth/grid breaks
  struct sfile_breaks
  {
    int p; // patch id
    int g; // grid # (0-based from bottom)
    int ib; // i index on grid g to start at
    int ie; // i index on grid g to end at
    int jb; // i index on grid g to start at
    int je; // i index on grid g to end at
    int kb; // k index on grid g to start at (1-based from top of grid)
    int ke; // k index on grid g to end at (1-based from top of grid)
    int hs; // horizontal sampling for this grid to this patch (>= 1)
    int vs; // vertical sampling for this grid to this patch (>= 1)
  };

  // Helper functions, just for internal use
  static void calculate_patches(EW& ew, vector<float_sw4>& mr_depth, 
      int horizontalSampling, vector<vector<sfile_breaks> >& breaks,
      vector<int>& patch_nk);

#ifdef USE_HDF5
  static void write_sfile_header(hid_t file_id, hid_t mpiprop_id,
      const float& h_coarse, const float (&lonlataz)[3], vector<int>& patch_nk);
  static void write_sfile_interfaces(hid_t file_id, hid_t mpiprop_id,
      EW& ew, vector<vector<sfile_breaks> >& breaks, 
      vector<float*>& z_bot, vector<float*>& z_top);
  static void write_sfile_materials(hid_t file_id, hid_t mpiprop_id, EW& ew,
      MaterialRfile& model, vector<vector<sfile_breaks> >& patch_breaks,
      vector<int>& patch_nk, vector<float*>& z_bot, vector<float*>& z_top);
  static void material_interpolate(vector<float*>& h5_array,
      float* zbot, float* ztop, hsize_t (&slice_dims)[3], 
      sfile_breaks& brk, Sarray* z, float zmin, float gridh,
      float& Cs_min, float& Cs_max,
      Sarray& rho, Sarray& mu, Sarray& lambda, Sarray& qp, Sarray& qs);
  static void patch_interface(float* z, hsize_t (&dims)[2], 
      bool top, vector<sfile_breaks>& pbrk, EW& ew);
  static void read_sfile_header(hid_t file_id, hid_t mpiprop_id, 
      float& h, float (&lonlataz)[3], vector<int>& patch_nk);
  static void read_sfile_interface_group(hid_t file_id, hid_t mpiprop_id, 
      bool topoOnly, vector<Sarray*>& intf);
  static void read_sfile_material_group(hid_t file_id, hid_t mpiprop_id, 
      int nghost, int nvars, vector<Sarray>& matl);
  static void calculate_grid_boundingbox(EW& ew, float_sw4 (&bb)[3][2]);
  static void calculate_interpolation_patch(vector<Sarray>& matl, 
      int nghost, float_sw4 (&bb)[3][2], float_sw4 x0, float_sw4 m_y0, 
      float hh, int nvars, vector<int>& patch_nk);
  static void get_patch_dims( sfile_breaks brk, int& ibeg, int& iend, int& jbeg, int& jend );
#endif // ifdef USE_HDF5

private:
  SfileHDF5(); // make it impossible to call default constructor
  SfileHDF5(const SfileHDF5 &in); // hide copy constructor
  void operator=(const SfileHDF5 &in); // hide assignment operator
};

#endif
