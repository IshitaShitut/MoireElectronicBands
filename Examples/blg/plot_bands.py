from pymelecutil.pymelecutil import plot_data

label=(r'$\Gamma $','M', 'K', r'$\Gamma $')
en_range = [-15,8.5]
plt = plot_data(dpi=500,en_range=en_range)
plt.band_structure(data_file='bands.dat.hdf5',
                   label=label, 
                   nbands=4, 
                   nkpt=300, 
                   kfile='k_points.dat',
                   closed_loop=True,
                   save=False)
