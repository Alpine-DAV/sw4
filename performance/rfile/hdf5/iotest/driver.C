#include <mpi.h>
#include "EW.h"
#include "Require.h"

// Simple driver to test hdf5 parallel
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int myRank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

  ASSERT(argc == 2);
  // string filename = "berkeley-rfile-2-sfile.in";
  string filename(argv[1]);

  vector<vector<Source*> > vecvec_source;
  vector<vector<TimeSeries*> > vecvec_series;
  // EW ew(filename, vecvec_source,  vecvec_series);
  EW* ew = NULL;
  string rfile_name = "berkeley.rfile";
  string rfile_dir = ".";
  int bufsize = 200000;
  // MaterialRfile mat_rfile(ew, rfile_name, rfile_dir, bufsize);

  string sfile_name = "berkeley.sfile";
  vector<double> vec_depths;
  // Hardwire for berkeley.rfile
  vec_depths.push_back(500);
  vec_depths.push_back(2000);
  // vec_depths.push_back(6387.5);
  MaterialSfile mat_sfile(ew, sfile_name, rfile_dir, false, true, 1, vec_depths);
  
  mat_sfile.write_sfile(rfile_dir, rfile_name);
  /*
  MaterialRfile mat_rfile;
  mat_rfile.setup();
  SfileHDF5 writer("berkeley.rfile", ".");
  writer.write_sfile(ew1, mat_rfile);

  cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  EW ew2;
  ew2.setup();
  writer.read_sfile_topo(ew2);
  MaterialSfile mat_sfile;
  writer.read_sfile_mat(ew2, mat_sfile);
  */

  MPI_Finalize();

  return 0;
}
