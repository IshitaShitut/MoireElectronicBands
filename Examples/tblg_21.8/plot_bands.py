from pymelecutil.pymelecutil import plot_data

label=(r'$\Gamma $','M', 'K', r'$\Gamma $')
en_range = [-13,7.5]
plt = plot_data(dpi=500,en_range=en_range)
plt.band_structure(data_file='bands_21.8_E0.hdf5',
                   label=label, 
                   nbands=28, 
                   nkpt=300, 
                   kfile='k_points.dat',
                   closed_loop=True,
                   save=False)
