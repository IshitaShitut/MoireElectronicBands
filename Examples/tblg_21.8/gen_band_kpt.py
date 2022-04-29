from pymelecutil.pymelecutil import moire_electron_utils
import numpy as np

gamma = np.array([0,0,0])
K = np.array([1/3.,1/3.,0])
M = np.array([0.5,0.0,0])
Kp = np.array([2/3,-1/3,0])

number_of_points_path_1 = 10
number_of_points_path_2 = 8
number_of_points_path_3 = 8
number_of_points_path_4 = 10

path = [gamma, number_of_points_path_1, K,
        number_of_points_path_2, M,
        number_of_points_path_3, Kp,
        number_of_points_path_4, gamma]

points = moire_electron_utils()
points.get_kpt_bands(path,output_file='k_points_test.dat')
