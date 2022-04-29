from pymelecutil.pymelecutil import plot_data

label=(r'$\Gamma $','K', 'M','K\'',r'$\Gamma $')
en_range = [-1.000,1.000]
plt = plot_data(dpi=600,en_range=en_range)
plt.band_structure(data_file='bands_grhbn_relaxed_v2.hdf5',
                   label=label,
                   nbands=1001, 
                   nkpt=36,
                   kfile='k_points_test.dat',
                   closed_loop=True,
                   save=False)

