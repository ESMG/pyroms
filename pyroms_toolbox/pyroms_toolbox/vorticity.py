import numpy as np

def vorticity(u, v, grd):
    """
    compute the relative vorticity
    """

    dx = grd.hgrid.dx
    dy = grd.hgrid.dy

    #dx, dy at psi point
    dx = 0.5 * (dx[1:,:] + dx[:-1,:])
    dy = 0.5 * (dy[:,1:] + dy[:,:-1])

    vorticity = np.zeros(grd.hgrid.mask_psi.shape)

    vorticity = 2 * (v[:,1:] - v[:,:-1]) / (dx[:,1:] + dx[:,:-1]) \
                - 2 * (u[1:,:] - u[:-1,:]) / (dy[1:,:] + dy[:-1,:])

    vorticity = np.ma.masked_where(grd.hgrid.mask_psi == 0, vorticity)

    return vorticity
