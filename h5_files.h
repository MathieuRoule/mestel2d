#include <valarray>
#include <string>
#include <hdf5.h>
#include "hdf5_hl.h"
using namespace std;

/*
 *
 *  The base class
 *
 */
class H5base {
protected:
  hid_t id;
public:
  H5base() {};
  ~H5base() {};
  hid_t operator()() { return id; }
  void write(real x, const string &name,
	     const string &title="");
  void write(valarray<real> &x, const string &name,
	     const string &title="");
  void write(valarray<real> &x, const int ndim, const hsize_t *dims,
       const string &name, const string &title="");
  herr_t set(const string &name, const string &value);
  herr_t set(const string &name, float value);
  herr_t set(const string &name, double value);
  herr_t set(const string &name, int value);
  real getnum(const string &name);
  std::string getstring(const string &name);
  std::valarray<real> readvec(const string &name);
};


/*
 *
 *  The file class
 *
 */
class H5file : public H5base {
  public:
  H5file(const string &fnam, const char mode='r',
          const string &title="Grommet file");
  ~H5file() { if(id!=0) H5Fclose(id); id=0; }
  herr_t close() { herr_t status=H5Fclose(id); id=0; return status; }
  herr_t flush() { return H5Fflush(id, H5F_SCOPE_LOCAL);  }
  // default copy constructor?
};

/*
 *
 *  The group class
 *
 */
class H5group : public H5base {
public:
  H5group() {id=0; }
  H5group(const H5group &in) { id=in.id;}
  H5group(H5base &loc, const string &groupname, const string &title="");
  ~H5group() {
    id=0;
  }
  herr_t close() {
    if(id!=0) return H5Gclose(id);
    id=0; return 0;
  }
  bool exists(H5base &base, const string &groupname);
};
