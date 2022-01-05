import h5py
import numpy as np

fname = 'bands_21.8_E0-offset_v20.hdf5'

data = h5py.File(fname, 'r')

for group in data.keys():
    ds_data = data[group]['evec_real']
    evec = ds_data[:]
    for i in range(4):
        print("/n")
        for j in range(len(evec[i,:])):
            print(evec[i,j])

