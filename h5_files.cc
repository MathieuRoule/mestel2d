#include "h5_files.h"

using namespace std;

herr_t H5base::set(const string &name, const string &value) {
  return H5LTset_attribute_string( id, ".", name.c_str(), value.c_str() );
}

herr_t H5base::set(const string &name, float value) {
  return H5LTset_attribute_float( id, ".", name.c_str(), &value, 1 );
}
herr_t H5base::set(const string &name, double value) {
  return H5LTset_attribute_double( id, ".", name.c_str(), &value, 1 );
}
herr_t H5base::set(const string &name, int value) {
  return H5LTset_attribute_int( id, ".", name.c_str(), &value, 1 );
}
std::string H5base::getstring(const string &name) {
  char data[999];
  H5LTget_attribute_string(id, ".", name.c_str(), data);
  return std::string(data);
}
real H5base::getnum(const string &name) {
  double data;
  H5LTget_attribute_double(id, ".", name.c_str(), &data);
  return data;
}

valarray<real> H5base::readvec(const string &name) {
  hid_t dataset_id = H5Dopen(id,name.c_str());
  hid_t dataspace_id = H5Dget_space(dataset_id);
  int n = (int) H5Sget_simple_extent_npoints(dataspace_id);
  real *tmpbuf = new real[n];
  H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL,
          H5S_ALL, H5P_DEFAULT, tmpbuf);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);
  valarray<real> ans(n);
  for(int i=0;i<n;i++) ans[i]=tmpbuf[i];
  delete[] tmpbuf;
  return ans;
}

void H5base::write(real x,
                    const string &name, const string &title) {

  hsize_t dims[]={1};
  real *tmpbuf = new real[1];
  tmpbuf[0] = x;
  hid_t dataspace_id = H5Screate_simple(1,dims,NULL);
  if(H5LTfind_dataset ( id, name.c_str() )) {
    H5Gunlink(id,name.c_str());
  }
  hid_t dataset_id = H5Dcreate(id, name.c_str(), H5T_NATIVE_FLOAT, 
                               dataspace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT,
			   H5S_ALL, H5S_ALL,H5P_DEFAULT, tmpbuf);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  delete[] tmpbuf;
}

void H5base::write(valarray<real> &arr,
                    const string &name, const string &title) {
  int n = arr.size();
  hsize_t dims[]={(hsize_t) n,};
  real *tmpbuf = new real[n];
  for(int i=0;i<n;i++) tmpbuf[i] = arr[i];
  hid_t dataspace_id = H5Screate_simple(1,dims,NULL);
  if(H5LTfind_dataset ( id, name.c_str() )) 
    H5Gunlink(id,name.c_str());
  hid_t dataset_id = H5Dcreate(id, name.c_str(), H5T_NATIVE_FLOAT, 
                               dataspace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT,
			   H5S_ALL, H5S_ALL,H5P_DEFAULT, tmpbuf);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  delete[] tmpbuf;
}


void H5base::write(valarray<real> &arr, const int ndim, const hsize_t *dims,
                    const string &name, const string &title) {
  int ntot = 1;
  for (int d=0;d<ndim;d++) ntot *= dims[d];
  real *tmpbuf = new real[ntot];
  for(int i=0;i<ntot;i++) tmpbuf[i] = arr[i];
  hid_t dataspace_id = H5Screate_simple(ndim,dims,NULL);
  if(H5LTfind_dataset ( id, name.c_str() )) 
    H5Gunlink(id,name.c_str());
  hid_t dataset_id = H5Dcreate(id, name.c_str(), H5T_NATIVE_FLOAT, 
                               dataspace_id, H5P_DEFAULT);

  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT,
			   H5S_ALL, H5S_ALL,H5P_DEFAULT, tmpbuf);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  delete[] tmpbuf;
}



H5file::H5file(const string &fnam, const char mode, const string &title) {
  if(mode=='r')  // open for reading only
    id = H5Fopen(fnam.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  else if(mode=='n') {// open new file
    id = H5Fcreate(fnam.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  } else if(mode=='w') {// open new file, truncating
    id = H5Fcreate(fnam.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }  else if(1||mode=='a') {
    H5E_auto_t err_fun;
    void *err_data;
    H5Eget_auto(&err_fun, &err_data);
    H5Eset_auto(NULL,NULL);
    id = H5Fcreate(fnam.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if(id<0) id = H5Fopen(fnam.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    H5Eset_auto(err_fun, err_data);
  }
}

bool H5group::exists(H5base &base, const string &groupname) {
  hid_t loc = base();
  H5E_auto_t err_fun;
  void *err_data;
  H5Eget_auto(&err_fun, &err_data);
  herr_t status = H5Eset_auto(NULL, NULL);
  status = H5Gget_objinfo (loc, groupname.c_str(), 0, NULL);
  bool found = (status==0) ? true : false;
  H5Eset_auto(err_fun, err_data);
  return found;
}


H5group::H5group(H5base &base, const string &groupname,
                   const string &title) {
  hid_t loc = base();
  if(exists(base, groupname.c_str()))
    id = H5Gopen(loc,groupname.c_str());
  else {
    id = H5Gcreate(loc, groupname.c_str(), H5P_DEFAULT);
  }
}
