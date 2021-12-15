import numpy as np

def k_pt_bands(points):
    """
        Generates a list of points along the path specified for plotting the 
        band structure

        Input: 
              points: list of k-points along which the bands are plotted
        
              The input style of points is - 
              [..., k(i), no of points between k(i) 
                                        and k(i+1), k(i+1), ...] 
              where k(n) is the n-th k-point, in fractional coordinates
                                            
        
        Output: 
                k_pt: list of the k-point path in fractional coordinates
    """
    k_pt, nodes = [], [0,]
    for i in range(0,len(points)-1,2):
        ki = points[i]
        kiplus1 = points[i+2]
        num = points[i+1]
        nodes.append(num+nodes[len(nodes)-1])
        for i in range(num):
            k_pt.append(ki+((kiplus1-ki)*i/num))
    k_pt.append(points[len(points)-1])
    return np.array(k_pt), nodes


def __main__():
    gamma = np.array([0,0,0])
    M = np.array([0.5,0.0,0])
    K = np.array([2/3.,1/3.,0])

    number_of_points_path_1 = 10
    number_of_points_path_2 = 8
    number_of_points_path_3 = 12

    path = [gamma, number_of_points_path_1, M,
            number_of_points_path_2, K,
            number_of_points_path_3, gamma]

    k_point, nodes = k_pt_bands(path)

    f = open('k_points.dat','w')
    for i in range(len(k_point)):
        f.write('%+.6f\t%+.6f\t%+.6f\n'%(k_point[i,0],k_point[i,1],k_point[i,2]))
    f.close()


__main__()
