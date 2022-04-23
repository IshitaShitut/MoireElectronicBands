from pymelecutil.pymelecutil import plot_data

label=(r'$\Gamma $','M', 'K', r'$\Gamma $')
en_range = [-0.250,0.250]
plt = plot_data(dpi=500,en_range=en_range)
#plt.band_structure(data_file='../../Examples/tblg_21.8/bands_21.8_E0.hdf5',label=label, 
#                        nbands=28, nkpt=300, kfile='../../Examples/tblg_21.8/k_points.dat',
#                        save=False)
plt.band_structure(data_file='../../Examples/tblg_1.61/bands_1.61_E0.hdf5',label=label, 
                   nbands=44, nkpt=31, kfile='../../Examples/tblg_1.61/k_points.dat',
                   closed_loop=False,save=False)
