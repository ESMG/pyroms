import numpy as np

def rx1(z_w,rmask):
    """
    function rx1 = rx1(z_w,rmask)
 
    This function computes the bathymetry slope from a SCRUM NetCDF file.

    On Input:
       z_w         layer depth.
       rmask       Land/Sea masking at RHO-points.
 
    On Output:
       rx1         Haney stiffness ratios.
    """

    N, Lp, Mp = z_w.shape
    L=Lp-1
    M=Mp-1

    #  Land/Sea mask on U-points.
    umask = np.zeros((L,Mp))
    for j in range(Mp):
        for i in range(1,Lp):
            umask[i-1,j] = rmask[i,j] * rmask[i-1,j]

    #  Land/Sea mask on V-points.
    vmask = np.zeros((Lp,M))
    for j in range(1,Mp):
        for i in range(Lp):
            vmask[i,j-1] = rmask[i,j] * rmask[i,j-1]

    #-------------------------------------------------------------------
    #  Compute R-factor.
    #-------------------------------------------------------------------

    zx = np.zeros((N,L,Mp))
    zy = np.zeros((N,Lp,M))

    for k in range(N):
        zx[k,:] = abs((z_w[k,1:,:] - z_w[k,:-1,:] + z_w[k-1,1:,:] - z_w[k-1,:-1,:]) / 
                      (z_w[k,1:,:] + z_w[k,:-1,:] - z_w[k-1,1:,:] - z_w[k-1,:-1,:]))
        zy[k,:] = abs((z_w[k,:,1:] - z_w[k,:,:-1] + z_w[k-1,:,1:] - z_w[k-1,:,:-1]) /
                      (z_w[k,:,1:] + z_w[k,:,:-1] - z_w[k-1,:,1:] - z_w[k-1,:,:-1]))
        zx[k,:] = zx[k,:] * umask
        zy[k,:] = zy[k,:] * vmask


    r = np.maximum(np.maximum(zx[:,:,:-1],zx[:,:,1:]), np.maximum(zy[:,:-1,:],zy[:,1:,:]))

    rx1 = np.amax(r, axis=0)

    rmin = rx1.min()
    rmax = rx1.max()
    ravg = rx1.mean()
    rmed = np.median(rx1)

    print('  ')
    print('Minimum r-value = ', rmin)
    print('Maximum r-value = ', rmax)
    print('Mean    r-value = ', ravg)
    print('Median  r-value = ', rmed)

    return rx1
