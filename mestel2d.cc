#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <ctime>
#include <cmath>
#include "2d_box.h"

const real G = 1.0;

void usage(char *progname) {
  cout << "Usage: " << progname << " [options] ic-file dump-file" << endl;
  cout << "with options: " << endl;
  cout << "  -xmax=26       potential box encloses |x_i|<xmax" << endl
       << "  -nx=120        with nx cells across" << endl
       << "  -dt=1e-3       integrate with timestep dt" << endl
       << "  -partdumpevery=200 dumping summary stat to dump-file every partdumpevery steps" << endl
       << "  -fulldumpevery=200 dumping a snapshot to dump-file every fulldumpevery steps" << endl
       << "  -nstep=1000000 for nstep steps" << endl
       << "  -nring=-1      with nring radial cells (if =-1, default value is 20*nx), only used if m2=true" << endl
       << "  -nphi=720      with nphi azimuthal cells, only used if m2=true" << endl
       << "  -kernel='plummer' gravity softening kernel (plummer/kuzmin), default is plummer" << endl
       << "  -eps=0.18      smmoothing of the interaction potential" << endl
       << "  -tstop=120     final time (if <0 uses nstep)" << endl
       << "  -fractive=1.0  active fraction" << endl
       << "  -m2            filter out all but m=2 harmonics" << endl
       << "  -basis='gate'  radial basis functions for coefficient computations (gate/gaussian/lognorm)" << endl
       << "  -rb=3.0        characteristic radius of basis functions" << endl
       << "  -sigmar=2.0    radial extension of basis functions" << endl
       << "  -lmax=2        maximal azimuthal number for the basis functions" << endl
       << "  -verbose=0     verbose"
       << endl
;
  exit(1);
}
int main(int argc, char *argv[]) {
  real vc2 = 1.0;                            // Circular velocity (squared): V0^2.
  real tstop = -1, dt=1e-3;                  // Final time (if <0 uses nstep) and time step.
  int nstep = 1000000;                       // Number of timesteps (used if tstop<0).
  int partdumpevery=10, fulldumpevery=200;   // Coefficients and full dumping rate.
  real xmax = 26.0;                          // Set the box size: (2*xmax)^2. Outside the box, only Mestel potential.
  int nx = 120;                              // Set the box cells size: (2*xmax/nx).
  string kernel = "plummer";                 // Gravity softening kernel, default: plummer softening.
  real eps = 0.18;                           // Smoothing of the interaction potential. Should be of order the grid size.
  real fractive = 1.0;                       // Active fraction (simple mass prefactor at the very early stage -- usefull to keep the same ICs).
  bool m2 = false;                           // m2=false -> Mestel potential + all m>0; m2=true -> Mestel potential + m=2.
  int nring = -1, nphi = 720;                // For azimuthal grid (if =-1, default value is 20*nx), only used if m2=true.
  string basis = "gate";                     // Basis functions used for coefficient computations.
  real rb = 3.0, sigmar = 2.0;               // Characteric radius and width for basis functions.
  int lmax = 2;                              // Maximal harmonic considered for basis functions.
  int verbose = 0;                           // Verbose level.
  
  int opt;
  for(opt=1; opt<argc; opt++) { 
    if (argv[opt][0] != '-') break;
    int ok = 0, p;
    char *option = new char[1+strlen(argv[opt])];
    strcpy(option,argv[opt]);
    for(p=1; option[p]!='=' && p<(int)strlen(option) ; p++) ;
    option[p] = 0;
    char *value = option+p+1;
    if (!strcmp(option, "-xmax"))   { xmax = atof(value); ok=1; }
    if (!strcmp(option, "-nx"))     { nx = atoi(value); ok=1; }
    if (!strcmp(option, "-dt"))     { dt = atof(value); ok=1; }
    if (!strcmp(option, "-partdumpevery")) { partdumpevery = atoi(value); ok=1; }
    if (!strcmp(option, "-fulldumpevery")) { fulldumpevery = atoi(value); ok=1; }
    if (!strcmp(option, "-nstep")) { nstep = atoi(value); ok=1; }
    if (!strcmp(option, "-kernel"))  { kernel=value; ok=1; }
    if (!strcmp(option, "-eps"))   { eps = atof(value); ok=1; }
    if (!strcmp(option, "-tstop")) { tstop = atof(value); ok=1; }
    if (!strcmp(option, "-fractive")) { fractive = atof(value); ok=1; }
    if (!strcmp(option, "-m2"))  { m2=true; ok=1; }
    if (!strcmp(option, "-nring")) { nring = atoi(value); ok=1; }
    if (!strcmp(option, "-nphi"))  { nphi = atoi(value); ok=1; }
    if (!strcmp(option, "-basis"))  { basis=value; ok=1; }
    if (!strcmp(option, "-rb")) { rb = atof(value); ok=1; }
    if (!strcmp(option, "-sigmar")) { sigmar = atof(value); ok=1; }
    if (!strcmp(option, "-lmax"))  { lmax = atoi(value); ok=1; }
    if (!strcmp(option, "-verbose"))  { verbose = atoi(value); ok=1; }
    if(!ok) {
      cerr << "I don't understand " << argv[opt] << "!" << endl;
      exit(1);
   }
    delete[] option;
  }
  if(argc-opt!=2) usage(argv[0]);
  string infile = string(argv[opt++]);
  string outfile = string(argv[opt++]);
  if(tstop>0) nstep = (int) (tstop/dt);
  if(nring<0) nring = nx*20; // choose ~20 rings per grid cell for m=2 restriction
  
  // Read ICs
  H5file fpin(infile,'r');
  Bodies ptle(fpin,"/");
  fpin.close();
  int N = ptle.n;
  for(int n=0;n<N;n++) ptle.m[n] *= fractive;
  
  // Box for potential
  Box box(nx,xmax,eps,G,kernel);
  
  // Create the output file
  H5file fout(outfile,'w');
  // Writing the run date
  std::chrono::time_point<std::chrono::system_clock> date = std::chrono::system_clock::now();
  std::time_t date_time = std::chrono::system_clock::to_time_t(date);
  fout.set("time",std::ctime(&date_time));
  // Dumping the parameters
  H5group params(fout,"params");
  params.set("xmax",xmax);
  params.set("nx",nx);
  params.set("dt",dt);
  params.set("partdumpevery",partdumpevery);
  params.set("fulldumpevery",fulldumpevery);
  params.set("nstep",nstep);
  params.set("kernel",kernel);
  params.set("eps",eps);
  params.set("fractive",fractive);
  if (m2) {
    params.set("m2",1);
    params.set("nring",nring);
    params.set("nphi",nphi);
  }
  params.close();

    // Create /full group in output file.  We dump the full position and velocities under this.
  H5group grpfull = H5group(fout,"full");
  grpfull.close();
  fout.flush();
  fout.close();

  // Sizes for the coefficients and times storing
  const int nfulldump = (fulldumpevery==0) ? 0 : (int) (1+nstep/fulldumpevery); // Number of full dump
  const int npartdump = (partdumpevery==0) ? 0 : (int) (1+nstep/partdumpevery); // Number of coefficient (or partial) dump
  // WARNING: assumption ncoef == 2*(lmax+1), has to be consistent
  const int ncoef = 2*(lmax+1);                    // Total number of coefficients
  // Times storing
  valarray<real> tfulldump(nfulldump);
  valarray<real> tpartdump(npartdump);
  // Coefficient storing
  valarray<real> instantcoef(ncoef);    // To compute the coefficients at a given time
  valarray<real> Acoefs(npartdump*ncoef);   // To store all the coefficients at all times
  // Unfortunately dumping (dynamically allocated) multi-dimensional array
  // in HDF5 files has no simple solution -> easiest solution is to handle 
  // index and use a long one-dimensional array.
  hsize_t Adims[] = {(hsize_t) npartdump, (hsize_t) ncoef};
  // Counters
  int fullcount = 0, partcount = 0;

  // Now the leapfrog steps
  real t=0.0;
  for(int istep=0; istep<=nstep; istep++) {

    if((npartdump!=0)&&(istep%partdumpevery==0)) {
      // Storing the time
      tpartdump[partcount] = t;
      // Computing the coefficients
      if (basis=="gaussian") instantcoef = ptle.compute_Ap_gaussian(rb,sigmar,lmax);
      else if (basis=="lognorm") instantcoef = ptle.compute_Ap_lognorm(rb,sigmar,lmax);
      else instantcoef = ptle.compute_Ap_gate(rb-sigmar,rb+sigmar,lmax); // default="gate" (rmin,rmax,lmax)
      // Storing them
      for (int l=0;l<=lmax;l++) {
        Acoefs[partcount*ncoef+2*l] = instantcoef[2*l];
        Acoefs[partcount*ncoef+2*l+1] = instantcoef[2*l+1];
      }
      // Output message
      if(verbose>1) {
        cout << "Coefficient computed: " << (partcount+1) << " over " << npartdump << endl;
      }
      // Counter update
      partcount += 1;
    }

    if((nfulldump!=0)&&(istep%fulldumpevery==0)) {
      // Open
      H5file fout(outfile,'a');
      H5group grpfull = H5group(fout,"full");
      // Storing the time
      tfulldump[fullcount] = t;
      // Dumping the full positions and velocities
      char out[50];
      snprintf(out,50,"snap%04d",istep/fulldumpevery); // Dump under /runs/snap0000, /runs/snap0001 etc
      if (istep==0) ptle.write_all(grpfull,t,out); // For ics also dump the masses
      else ptle.write_xv(grpfull,t,out); // Do not dump the masses
      grpfull.close();
      fout.set("latest",out);
      fout.flush();
      fout.close();
      // Output message
      if(verbose>0) {
        cout << "Snapshot written to node /full/" << out << " in " << outfile << endl;
      }
      // Counter update
      fullcount += 1;
    }

    // Drift all ptles half a timestep
    ptle.xy += ptle.vxvy*dt/2;
    t += dt/2;

    // Kick
    box.zero();
    if(m2) box.bodies2density_m2(ptle,nring,nphi);
    else box.bodies2density(ptle);
    box.density2pot();
    valarray<real> accels = box.pot2accels(ptle);
    for(int i=0;i<ptle.n;i++) {
      real x = ptle.xy[2*i], y = ptle.xy[2*i+1];
      real R2 = x*x+y*y;
      real ax = accels[2*i], ay = accels[2*i+1];
      if(!m2) {// Subtract m=0 component
        real Rdota = x*ax+y*ay;
        ax -= Rdota*x/R2;
        ay -= Rdota*y/R2;
      }
      // Contribution from monopole
      // circular speed R*dPhi/dR = vc^2 = const
      ax -= vc2*x/R2;
      ay -= vc2*y/R2;
      accels[2*i] = ax;
      accels[2*i+1] = ay;
    }
    ptle.vxvy += accels*dt;
    
    // Drift another half timestep
    ptle.xy += ptle.vxvy*dt/2;
    t += dt/2;

  }
  // Dumping the full dump times in a single array
  if (nfulldump!=0){
    // Open
    H5file fout(outfile,'a');
    H5group grpfull = H5group(fout,"full");
    // Times
    grpfull.write(tfulldump,"tabt");
    // Close
    grpfull.close();
    fout.flush();
    fout.close();
    // Output message
    if(verbose>0) {
      cout << "Snapshots time list written in " << outfile << endl;
    }
  }
  // Dumping the coefficients values and times
  // Create /coefdump group in output file.  We dump the summary statistic (lharmonic coefficients) under this.
  if (npartdump!=0){
    // Open
    H5file fout(outfile,'a');
    H5group grpcoef = H5group(fout,"coef");
    // Parameters
    grpcoef.set("basis",basis);
    grpcoef.set("rb",rb);
    grpcoef.set("sigmar",sigmar);
    grpcoef.set("lmax",lmax);
    // Times
    grpcoef.write(tpartdump,"tabt");
    // Values
    grpcoef.write(Acoefs,2,Adims,"tabcoef");
    // Close
    grpcoef.close();
    fout.flush();
    fout.close();
    // Output message
    if(verbose>0) {
      cout << "Coefficients/summary written in " << outfile << endl;
    }
  }
}


