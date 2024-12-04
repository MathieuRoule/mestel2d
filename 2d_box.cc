#include "2d_box.h"
#include <iostream>
const real pi = M_PI;

Box::Box(int n, real xmax, real eps, real G, string kernel)
  : G(G), nx(n), ny(n), eps(eps), kernel(kernel),
    xmin(-xmax), ymin(-xmax),
    dx(2*xmax/n), dy(2*xmax/n),
    odx(1.0/dx), ody(1.0/dy),
    xmin1(xmin+dx), ymin1(ymin+dy),
    xmin2(xmin+2*dx), ymin2(ymin+2*dy),
    xmax2(-xmin2), ymax2(-ymin2)
{
// Get ready for poisson
  // FFTW requires strange array layout: see
  // http://fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
  fftw_mesh = (double *)
    fftw_malloc((size_t) nx*4*ny*sizeof(double));
  fftw_ftkernel = (fftw_complex*)
    fftw_malloc(sizeof(fftw_complex)*nx*2*(ny+1));
  fftw_ftmesh = (fftw_complex*)
    fftw_malloc(sizeof(fftw_complex)*nx*2*(ny+1));
  // Build the FFTW "plans" for carrying out the transforms
  pk = fftw_plan_dft_r2c_2d(nx*2, ny*2,
                            fftw_mesh, fftw_ftmesh, FFTW_MEASURE);
  pkinv = fftw_plan_dft_c2r_2d(nx*2, ny*2,
                               fftw_ftmesh, fftw_mesh, FFTW_MEASURE);

  // Calculate the FT of our Poisson kernel.
  // First set up the kernel with the usual zero padding
  double scale = G/(nx*ny*4);
  for(int i=0;i<nx*ny*4;i++) fftw_mesh[i]=0.0;
  for(int i=-nx+1;i<nx;i++)
    for(int j=-ny+1;j<ny;j++) {
      double usq=(i*i*dx*dx+j*j*dy*dy)/(eps*eps);
      int i0=(i+2*nx)%(2*nx), j0=(j+2*ny)%(2*ny);
      int ij=i0*2*ny+j0;
      if (kernel == "kuzmin"){
        fftw_mesh[ij]=-scale*(3.0+2.5*usq+usq*usq)/(eps*pow(sqrt(1.0+usq),5));// Kuzmin softening
      } else {
        // John's version
        // double r=sqrt(eps*eps+i*i*dx*dx+j*j*dy*dy); // softened radius
        // fftw_mesh[ij]=-scale/r;
        // New version
        fftw_mesh[ij]=-scale/(eps*sqrt(1.0+usq)); // Plummer softening
      }
    }
  // QUESTION !! It does not seem to be zero-padded at all ! (It's padded with U expected values)
  // The (zero-padded) kernel is set up.  Calculate its FT.
  fftw_execute(pk);
  // Now copy fftw_ftmesh to fftw_ftkernel
  for(int i=0;i<nx*2*(ny+1);i++) {
    fftw_ftkernel[i][0]=fftw_ftmesh[i][0];
    fftw_ftkernel[i][1]=fftw_ftmesh[i][1];
  }
}

Box::~Box() {
  fftw_free(fftw_mesh);
  fftw_free(fftw_ftmesh);
  fftw_free(fftw_ftkernel);
}

void Box::zero() {
  for(int i=0;i<nx*4*ny;i++) fftw_mesh[i] = 0.0;
}

void Box::density2pot() {
  // transform density (in fftw_mesh) -> FT(density) (in fftw_ftmesh)
  fftw_execute(pk);

  // Multiply by FT poisson kernel
  fftw_complex *msh = fftw_ftmesh, *krnl = fftw_ftkernel;
  for(int i=0;i<2*nx;i++)
    for(int j=0;j<ny+1;j++) {
      int ij=i*(ny+1)+j;
      double tmp=msh[ij][0];
      // Re(z*z') = Re(z)*Re(z') - Im(z)*Im(z')
      msh[ij][0]=krnl[ij][0]*tmp-krnl[ij][1]*msh[ij][1];
      // Im(z*z') = Re(z)*Im(z') + Im(z)*Re(z')
      msh[ij][1]=krnl[ij][0]*msh[ij][1]+krnl[ij][1]*tmp;
    }

  // Transform back FT(potential) (in fftw_ftmesh) -> potential (in fftw_mesh)
  fftw_execute(pkinv);
}

void Box::bodies2density(const Bodies &ptle) {
  for(int n=0;n<ptle.n;n++) {
    real x = ptle.xy[2*n];
    real y = ptle.xy[2*n+1];
    real m = 0.5*ptle.m[n]; // 0.5 because we reflect
    const int DX=2*ny, DY=1;
    if(fabs(x)<xmax2 && fabs(y)<ymax2) {
      int ix = int(odx*(x-xmin));
      int iy = int(ody*(y-ymin));
      int ixy = ix*2*ny+iy;
      real fx = odx*(x-xmin)-ix;
      real fy = ody*(y-ymin)-iy;
      fftw_mesh[ixy]    += (1-fx)*(1-fy)*m;
      fftw_mesh[ixy+DX] += fx*(1-fy)*m;
      fftw_mesh[ixy+DY]    += (1-fx)*fy*m;
      fftw_mesh[ixy+DX+DY] += fx*fy*m;
      // Reflect
      x=-x; y=-y;
      ix = int(odx*(x-xmin));
      iy = int(ody*(y-ymin));
      ixy = ix*2*ny+iy;
      fx = odx*(x-xmin)-ix;
      fy = ody*(y-ymin)-iy;
      fftw_mesh[ixy]    += (1-fx)*(1-fy)*m;
      fftw_mesh[ixy+DX] += fx*(1-fy)*m;
      fftw_mesh[ixy+DY]    += (1-fx)*fy*m;
      fftw_mesh[ixy+DX+DY] += fx*fy*m;
    }
  }
}

valarray<real> Box::pot2accels(const Bodies &ptle) {
  const int DX=2*ny, DY=1;
  valarray<real> ans(2*ptle.n);
  ans = 0.0;
  for(int n=0;n<ptle.n;n++) {
    real x = ptle.xy[2*n];
    real y = ptle.xy[2*n+1];
    real m = ptle.m[n];
    if(fabs(x)<xmax2 && fabs(y)<ymax2) {
      int ix = int(odx*(x-xmin));
      int iy = int(ody*(y-ymin));
      int ixy = ix*2*ny+iy;
      real fx = odx*(x-xmin)-ix;
      real fy = ody*(y-ymin)-iy;
      double *mxy = fftw_mesh+ixy;
      ans[2*n+0] = -0.5*odx*
	( (1-fx)*(1-fy)*(mxy[+DX]-mxy[-DX])
	  +(1-fx)*fy*(mxy[+DX+DY]-mxy[-DX+DY])
	  +fx*(1-fy)*(mxy[+2*DX]-mxy[0*DX])
	  +fx*fy*(mxy[+2*DX+DY]-mxy[+DY]) );
      ans[2*n+1] = -0.5*ody*
	( (1-fx)*(1-fy)*(mxy[+DY]-mxy[-DY])
	  +(1-fx)*fy*(mxy[+2*DY]-mxy[0*DY])
	  +fx*(1-fy)*(mxy[+DX+DY]-mxy[+DX-DY])
	  +fx*fy*(mxy[+DX+2*DY]-mxy[+DX]) );
    }
  }
  return ans;
}

void Box::bodies2density_m2(const Bodies &ptle, int nring, int nphi) {
  // Fourier analysis: calculate rho_k in shells whose
  // inner boundaries are given by Rmin+i*dR;
  real Rmin = 0.0, Rmax = -xmin;
  real dR = (Rmax-Rmin)/nring;
  real soft2 = 1e-6;
  valarray<real> rhok_cos(nring), rhok_sin(nring);
  rhok_cos = 0.0; rhok_sin = 0.0;
  for(int n=0;n<ptle.n;n++) {
    real x = ptle.xy[2*n], y = ptle.xy[2*n+1];
    real R = sqrt(x*x+y*y+soft2);
    real cosphi = x/R, sinphi = y/R;
    real cos2phi = cosphi*cosphi-sinphi*sinphi;
    real sin2phi = 2*sinphi*cosphi;
    int ndx = ((int) ((R-Rmin)/dR));
    if(ndx>=0 && ndx<nring) {
      rhok_cos[ndx] += ptle.m[n]*cos2phi;
      rhok_sin[ndx] += ptle.m[n]*sin2phi;
    }
  }
  // pre-calculate trig constants
  valarray<real> cosphi(nphi), sinphi(nphi), cos2phi(nphi), sin2phi(nphi);
  for(int iphi=0;iphi<nphi;iphi++) {
    double phi = iphi*2*pi/nphi;
    cosphi[iphi] = cos(phi);
    sinphi[iphi] = sin(phi);
    cos2phi[iphi] = cos(2*phi);
    sin2phi[iphi] = sin(2*phi);
  }
  // now assign Fourier-reconstructed mass to the mesh, ring by ring
  for(int iring=0;iring<nring;iring++){
    real R = Rmin+iring*dR;
    real rhok_cos_here = 2*pi*rhok_cos[iring]/(nphi*pi);
    real rhok_sin_here = 2*pi*rhok_sin[iring]/(nphi*pi);
    Bodies ring(nphi);
    for(int iphi=0;iphi<nphi;iphi++) {
      ring.xy[2*iphi+0] = R*cosphi[iphi];
      ring.xy[2*iphi+1] = R*sinphi[iphi];
      ring.m[iphi] = rhok_cos_here*cos2phi[iphi]
	+rhok_sin_here*sin2phi[iphi];
    }
    bodies2density(ring);
  }
}
