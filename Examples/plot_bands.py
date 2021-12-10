import h5py
import matplotlib.pyplot as plt
import numpy as np

fname = 'results.h5'
natoms = 5044
no_k_pts = 5
index_ = [i for i in range(no_k_pts)]

data = h5py.File(fname, 'r')

eigvals = np.empty((no_k_pts, natoms))

counter = 0
for group in data.keys():
    for dset in data[group].keys():
        ds_data = data[group][dset]
        eigvals[counter] = ds_data[:]
    counter+=1
print(counter)
x = [i for i in range(no_k_pts)]
for i in range(natoms):
    plt.plot(x, eigvals[:,i],c='k')
plt.show()
