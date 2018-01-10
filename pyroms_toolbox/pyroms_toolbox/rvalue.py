import numpy as np

def rvalue(h):
    """
    function rv = rvalue(h)

    This function compute ROMS stiffness parameter 
    (ripped from John Wilkin rvalue.m)

    On Input:
       h           bathymetry at RHO-points.
 
    On Output:
       rv          ROMS stiffness parameter.
    """
    #check that h is 2D
    if (len(h.squeeze().shape)!=2):
        raise ValueError('h must be two dimensions')

    #check whether h contains any NaNs
    if np.isnan(h).any(): raise Warning('the height array contains NaNs')

    #compute  diff(h)/2*mean(h) at each velocity grid point
    dhdx_u = np.diff(h, axis=1)
    dhdy_v = np.diff(h, axis=0)
    th_u = 2 * 0.5*(h[:,1:] + h[:,:-1])
    th_v = 2 * 0.5*(h[1:,:] + h[:-1,:])
    r_u = abs(dhdx_u / th_u)
    r_v = abs(dhdy_v / th_v)

    #for each rho point, find the maximum rvalue at the 4 
    #surrounding u,v points
    r_u = np.maximum(r_u[:,1:],r_u[:,:-1])
    r_v = np.maximum(r_v[1:,:],r_v[:-1,:])

    # pad rows and columns to give a result the same size as the input h
    r_u = np.c_[r_u[:,0], r_u, r_u[:,-1]]
    r_v = np.r_[np.tile(r_v[0,:], (1, 1)), r_v, np.tile(r_v[-1,:], (1, 1))]

    rv = np.maximum(r_u,r_v)

    return rv
