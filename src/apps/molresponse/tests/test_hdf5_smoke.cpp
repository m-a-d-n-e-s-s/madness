// test_hdf5_smoke.cpp — P1 build-wiring validation for the io-hdf5 thread.
// Pure HDF5 C API: create a file, write a dataset, read it back, verify bitwise.
// Deliberately NO MADNESS function code — that is P2 (docs/30_io_hdf5_survey_and_plan.md).
#include <hdf5.h>
#include <cstdio>
#include <vector>

int main() {
  const char* fname = "test_hdf5_smoke.h5";
  const hsize_t n = 6;
  const std::vector<double> wbuf = {1.0, 2.5, -3.0, 4.25, 5.0, 6.5};

  // --- write ---
  hid_t file  = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) { std::fprintf(stderr, "H5Fcreate failed\n"); return 1; }
  hid_t space = H5Screate_simple(1, &n, nullptr);
  hid_t dset  = H5Dcreate2(file, "data", H5T_NATIVE_DOUBLE, space,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wbuf.data());
  H5Dclose(dset); H5Sclose(space); H5Fclose(file);

  // --- read back ---
  std::vector<double> rbuf(n, 0.0);
  hid_t file2 = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset2 = H5Dopen2(file2, "data", H5P_DEFAULT);
  H5Dread(dset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf.data());
  H5Dclose(dset2); H5Fclose(file2);

  // --- verify bitwise ---
  for (hsize_t i = 0; i < n; ++i) {
    if (rbuf[i] != wbuf[i]) {
      std::fprintf(stderr, "MISMATCH at %llu: wrote %g read %g\n",
                   (unsigned long long)i, wbuf[i], rbuf[i]);
      return 1;
    }
  }
  unsigned maj, min, rel;
  H5get_libversion(&maj, &min, &rel);
  std::printf("test_hdf5_smoke: PASS (HDF5 %u.%u.%u, round-tripped %llu doubles)\n",
              maj, min, rel, (unsigned long long)n);
  return 0;
}
