#include "2d_bodies.h"
#include <cmath>

Bodies::Bodies(int N) : m(N), xy(2*N), vxvy(2*N), n(N) {}
Bodies::~Bodies() {}

Bodies::Bodies(H5base &base, const string &path) :
  xy(base.readvec(path+"/xy")),
  vxvy(base.readvec(path+"/vxvy")),
  m(base.readvec(path+"/m")),
  n(m.size()) {}

// write out in HDF5 format, with extra stuff to make PyTables happy
void Bodies::write_all(H5base &base, const real t, const string &path) {
  H5group grp(base, path);
  grp.write(t, "t", "time");
  grp.write(m, "m", "masses");
  grp.write(xy, "xy", "positions");
  grp.write(vxvy, "vxvy", "velocities");
  grp.close();
}
// write out in HDF5 format, with extra stuff to make PyTables happy
void Bodies::write_xv(H5base &base, const real t, const string &path) {
  H5group grp(base, path);
  grp.write(t,"t","time");
  grp.write(xy, "xy", "positions");
  grp.write(vxvy, "vxvy", "velocities");
  grp.close();
}

// Decomposition coefficient using gate function as basis element
valarray<real> Bodies::compute_Ap_gate(const real rmin, const real rmax, const int lmax) {

  // Allocate result array
  // WARNING: assumption ncoef == 2*(lmax+1), has to be consistent
  valarray<real> ans(2*(lmax+1)); //
  // Initialize the coefficients values: = is equivalent to fill with valarray
  // (Does not fill beyond the size, see dynamical allocation and valarray)
  ans = 0.0;
  real mass,x,y,r,phi;
  // Loop over particles
  for (int i=0;i<n;i++) {
    mass = m[i];
    x = xy[2*i], y = xy[2*i+1];
    r = sqrt(x*x+y*y);
    // If the particle is between rmin and rmax
    if (rmin<=r<=rmax) {
        phi = atan2(y,x);
        // Loop over harmonic number from 0 to lmax
        for (int l=0;l<=lmax;l++){
          ans[2*l] += mass*cos(l*phi);
          ans[2*l+1] += mass*sin(l*phi);
        }
    }
  }
  return ans;
}

// Decomposition coefficient using gaussian function as basis element
valarray<real> Bodies::compute_Ap_gaussian(const real rb, const real sigmar, const int lmax) {

  // Allocate result array
  // WARNING: assumption ncoef == 2*(lmax+1), has to be consistent
  valarray<real> ans(2*(lmax+1)); //
  // Initialize the coefficients values: = is equivalent to fill with valarray
  // (Does not fill beyond the size, see dynamical allocation and valarray)
  ans = 0.0;
  real mass,x,y,deltar,phi;
  real psir;
  // Loop over particles
  for (int i=0;i<n;i++) {
    mass = m[i];
    x = xy[2*i], y = xy[2*i+1];
    deltar = sqrt(x*x+y*y)-rb;
    psir = exp(-deltar*deltar/(2.0*sigmar*sigmar)); // Not normalized (hopefully unimportant)
    phi = atan2(y,x);
    // Loop over harmonic number from 0 to lmax
    for (int l=0;l<=lmax;l++){
      ans[2*l] += mass*psir*cos(l*phi);
      ans[2*l+1] += mass*psir*sin(l*phi);
    }
  }
  return ans;
}

// Decomposition coefficient using log-normal distribution as basis element
valarray<real> Bodies::compute_Ap_lognorm(const real rb, const real sigmar, const int lmax) {

  // Allocate result array
  // WARNING: assumption ncoef == 2*(lmax+1), has to be consistent
  valarray<real> ans(2*(lmax+1)); //
  // Initialize the coefficients values: = is equivalent to fill with valarray
  // (Does not fill beyond the size, see dynamical allocation and valarray)
  ans = 0.0;

  real mass,x,y,deltar,phi;
  real psir;
  real softln = 1e-6; // Softening for the log
  // Loop over particles
  for (int i=0;i<n;i++) {
    mass = m[i];
    x = xy[2*i], y = xy[2*i+1];
    deltar = 0.5*log(x*x+y*y+softln*softln)-rb; // log(r) - mu where r is softened as sqrt(x^2+y^2+eps^2)
    psir = exp(-deltar*deltar/(2.0*sigmar*sigmar)); // Not normalized (hopefully unimportant)
    phi = atan2(y,x);
    // Loop over harmonic number from 0 to lmax
    for (int l=0;l<=lmax;l++){
      ans[2*l] += mass*psir*cos(l*phi);
      ans[2*l+1] += mass*psir*sin(l*phi);
    }
  }
  return ans;
}