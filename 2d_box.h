#include "2d_bodies.h"
#include <fftw3.h>

class Box {
 public:
  Box(int n, real xmax, real eps, real G, string kernel);
  ~Box();
  void zero();
  void bodies2density(const Bodies &ptle);
  void bodies2density_m2(const Bodies &ptle, int nring, int nphi);
  void density2pot();
  valarray<real> pot2accels(const Bodies &ptle);

 private:
  const real G;
  const int nx, ny;
  const real xmin, ymin;
  const real eps;
  const real dx, dy, odx, ody;
  const real xmin1, ymin1;
  const real xmin2, ymin2;
  const real xmax2, ymax2;
  const string kernel;

  double *fftw_mesh;
  fftw_complex *fftw_ftmesh, *fftw_ftkernel;
  fftw_plan pk, pkinv;

};
