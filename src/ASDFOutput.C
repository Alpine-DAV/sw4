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
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>

#include "ASDFOutput.h"
#include "EW.h"

#ifdef USE_ASDF
#include "ASDF_init.h"
#include "ASDF_write.h"
#endif

using namespace std;

ASDFOutput::ASDFOutput( EW* a_ew, std::string fileName):
  m_ew(a_ew),
  m_fileName(fileName)
{
}

void ASDFOutput::writeASDF( vector<TimeSeries*>& a_TimeSeries)
{
#ifdef USE_ASDF
  cout << "Writing ASDF file to: " << m_fileName << endl;

  ASDF_initialize_hdf5();

  stringstream filePrefix;
  string path = m_ew->getPath();
//building the file name...
  if( path != "." )
    filePrefix << path;
  filePrefix << m_fileName;

  hid_t file_id = ASDF_create_new_file(
      filePrefix.str().c_str(), MPI_COMM_WORLD);
  if (file_id < 0)
  {
    cout << "Could not open ASDF hdf5 file: " << filePrefix << "!!!" << endl;
    MPI_Abort(MPI_COMM_WORLD, file_id);
  }

  // Generic header and waveform group
  ASDF_write_string_attribute(file_id, "file_format", "ASDF");
  ASDF_write_string_attribute(file_id, "file_version", "0.0.1b");
  hid_t waveforms_grp = ASDF_create_waveforms_group(file_id);

  // Only support 3 xyz or vel components output for now
  const int num_waveforms=3;
  // Allocate char* buffers
  const int buff_size=100;
  char **waveform_names = (char **) malloc(num_waveforms*sizeof(char *));
  char *event_name = (char *) malloc(buff_size*sizeof(char));
  sprintf(event_name, "UNKNOWN");
  for (int i=0; i < num_waveforms; ++i)
    waveform_names[i] = (char *) malloc(buff_size*sizeof(char));

  // Allocate a temp that can hold our biggest time series data
  int max_size = 0;
  for (int ts=0; ts < a_TimeSeries.size(); ts++)
    max_size = max(max_size, a_TimeSeries[ts]->getNsteps());
  float *data = (float *) malloc(max_size*sizeof(float));

  // For each time series, create the station group
  for (int ts=0; ts < a_TimeSeries.size(); ts++)
  {
    // Skip if it's not our point
    bool myPoint = a_TimeSeries[ts]->myPoint();
    if (!myPoint)
      continue;

    string station_name = a_TimeSeries[ts]->getStationName();
    TimeSeries::receiverMode mode = a_TimeSeries[ts]->getMode();
    if (!(mode == TimeSeries::Displacement) && !(mode == TimeSeries::Velocity))
    {
      cout << "ASDF hdf5 outputs only displacement or velocity ... skipping: " 
        << station_name << endl;
      continue;
    }

    hid_t station_grp = 
      ASDF_create_stations_group(waveforms_grp, station_name.c_str());
    if (station_grp < 0)
    {
      cout << "Could not create ASDF hdf5 station group: " << station_name
        << " ... skipping." << endl;
      continue;
    }

    hid_t data_id[num_waveforms];
    int nsamples=a_TimeSeries[ts]->getNsteps();
    int start_time=a_TimeSeries[ts]->get_shift(); // TODO - right value?
    double sampling_rate=m_ew->getTimeStep(); // TODO - right value?
    char disp_names[3][2] = {"x","y","z"}; // 2 for ?\0
    char vel_names[3][3] = {"ux","uy","uz"}; // 3 for ??\0
    for (int i=0; i < num_waveforms; ++i)
    {
      if (mode == TimeSeries::Displacement)
        sprintf(waveform_names[i], "displacement_%s", disp_names[i]);
      else
        sprintf(waveform_names[i], "velocity_%s", vel_names[i]);
    }

    ASDF_define_waveforms(station_grp, num_waveforms, nsamples, start_time, 
      sampling_rate, event_name, waveform_names, data_id);

    assert(nsamples);
    int chunk_size = 10; // TODO - what's best write chunk?
    int num_chunk = (nsamples - 1) / chunk_size + 1;

    float_sw4** ts_data = a_TimeSeries[ts]->getRecordingArray();
    for (int c = 0; c < num_chunk; ++c)
    {
      int offset = chunk_size * c;
      int samples = (offset+chunk_size < nsamples) ? chunk_size :
        (nsamples % chunk_size);

      for (int i = 0; i < num_waveforms; ++i)
      {
        for (int j=0; j < samples; j++)
          data[j] = ts_data[i][j+offset];
        ASDF_write_partial_waveform(data_id[i], data, offset, samples);
      }
    }

    for (int i = 0; i < num_waveforms; ++i)
      H5Dclose(data_id[i]);

    ASDF_close_group(station_grp);
  }

  // Free mallocs
  for (int i = 0; i < num_waveforms; ++i)
    free(waveform_names[i]);
  free(data);
  free(waveform_names);
  free(event_name);

  // Close everything
  ASDF_close_group(waveforms_grp);
  H5Fclose(file_id);

  ASDF_finalize_hdf5();
#endif
}
