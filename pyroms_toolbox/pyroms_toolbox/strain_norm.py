import numpy as np

def strain_norm(u, v, grd):

    dx = grd.hgrid.dx
    dy = grd.hgrid.dy

    #u, v, dx, dy at psi point
    u = 0.5 * (u[1:,:] + u[:-1,:])
    v = 0.5 * (v[:,1:] + v[:,:-1])
    dx = 0.5 * (dx[1:,:] + dx[:-1,:])
    dx = 0.5 * (dx[:,1:] + dx[:,:-1])
    dy = 0.5 * (dy[:,1:] + dy[:,:-1])
    dy = 0.5 * (dy[1:,:] + dy[:-1,:])

    dudx = np.zeros(grd.hgrid.mask_psi.shape)
    dudy = np.zeros(grd.hgrid.mask_psi.shape)
    dvdx = np.zeros(grd.hgrid.mask_psi.shape)
    dvdy = np.zeros(grd.hgrid.mask_psi.shape)


    dudx[:,0] = (u[:,1] - u[:,0]) / dx[:,0]
    dudx[:,1:-1] = (u[:,2:] - u[:,:-2]) / (2 * dx[:,1:-1])
    dudx[:,-1] = (u[:,-1] - u[:,-2]) / dx[:,-1]

    dudy[0,:] = (u[1,:] - u[0,:]) / dy[0,:]
    dudy[1:-1,:] = (u[2:,:] - u[:-2,:]) / (2 * dy[1:-1,:])
    dudy[-1,:] = (u[-1,:] - u[-2,:]) / dy[-1,:]

    dvdx[:,0] = (v[:,1] - v[:,0]) / dx[:,0]
    dvdx[:,1:-1] = (v[:,2:] - v[:,:-2]) / (2 * dx[:,1:-1])
    dvdx[:,-1] = (v[:,-1] - v[:,-2]) / dx[:,-1]

    dvdy[0,:] = (v[1,:] - v[0,:]) / dy[0,:]
    dvdy[1:-1,:] = (v[2:,:] - v[:-2,:]) / (2 * dy[1:-1,:])
    dvdy[-1,:] = (v[-1,:] - v[-2,:]) / dy[-1,:]

    strain_norm = np.sqrt(dudx**2 + 0.25*dudy**2 + 0.25*dvdx**2 + dvdy**2)

    spval = -1e37
    idx = np.where(grd.hgrid.mask_psi == 0)
    strain_norm[idx] = spval
    idx = np.where(strain_norm > 1000)
    strain_norm[idx] = spval
    strain_norm = np.ma.masked_values(strain_norm, spval)

    return strain_norm
