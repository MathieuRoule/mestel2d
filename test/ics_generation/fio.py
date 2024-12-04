import numpy as np
import h5py

def load_ptle(fnam,where=''):
    fp = h5py.File(fnam,'r')
    m = np.array(fp[f'{where}/m'])
    x = np.array(fp[f'{where}/xy'])
    v = np.array(fp[f'{where}/vxvy'])
    fp.close()
    return m, x, v

def dump_ptle(m,x,v,fnam,where='',mode='w'):
    fp = h5py.File(fnam,mode)
    try: del fp[f'{where}/m']
    except KeyError: pass
    try: del fp[f'{where}/xy']
    except KeyError: pass
    try: del fp[f'{where}/vxvy']
    except KeyError: pass
    fp.create_dataset(f'{where}/m',data=m)
    fp.create_dataset(f'{where}/xy',data=x)
    fp.create_dataset(f'{where}/vxvy',data=v)
    fp.close()
    
