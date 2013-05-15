import numpy as np

def laplacian(field, dx, dy):
#   dx = grd.hgrid.dx
#   dy = grd.hgrid.dy
    rows, cols = field.shape

    mr, mc = field.shape
    laplacian = np.zeros ((mr, mc))

    ## x direction
    ## left and right boundary
    laplacian[:, 0] = (field[:, 0] - 2 * field[:, 1] + field[:, 2]) / (dx[:,0] * dx[:,1])
    laplacian[:, mc-1] = (field[:, mc - 3] - 2 * field[:, mc - 2] + field[:, mc-1]) \
        / (dx[:,mc - 3] * dx[:,mc - 2])

    ## interior points
    tmp1 = laplacian[:, 1:mc - 1]
    tmp2 = (field[:, 2:mc] - 2 * field[:, 1:mc - 1] + field[:, 0:mc - 2])
    tmp3 = dx[:,0:mc -2] * dx[:,1:mc - 1] 
    laplacian[:, 1:mc - 1] = tmp1 + tmp2 / tmp3
          
    ## y direction
    ## top and bottom boundary
    laplacian[0, :] = laplacian[0,:]  + \
        (field[0, :] - 2 * field[1, :] + field[2, :] ) / (dy[0,:] * dy[1,:])
      
    laplacian[mr-1, :] = laplacian[mr-1, :] \
        + (field[mr-3,:] - 2 * field[mr-2, :] + field[mr-1, :]) \
        / (dy[mr-3,:] * dy[mr-2,:]) 
    
    ## interior points
    tmp1 = laplacian[1:mr-1, :]
    tmp2 = (field[2:mr, :] - 2 * field[1:mr - 1, :] + field[0:mr-2, :])
    tmp3 = dy[0:mr-2,:] * dy[1:mr-1,:]
    laplacian[1:mr-1, :] = tmp1 + tmp2 / tmp3
    
    return laplacian

