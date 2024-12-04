#include "h5_files.h"

class Bodies {
 public:
  valarray<real> m, xy, vxvy;
  const int n;
  Bodies(int);
  Bodies(H5base &, const string &path);
  virtual ~Bodies();
  void write_all(H5base &, const real t, const string &path=".");
  void write_xv(H5base &, const real t, const string &path=".");
  valarray<real> compute_Ap_gate(const real rmin, const real rmax, const int lmax);
  valarray<real> compute_Ap_gaussian(const real rb, const real sigmar, const int lmax);
  valarray<real> compute_Ap_lognorm(const real rb, const real sigmar, const int lmax);
};
