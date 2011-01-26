import numpy as np

def rx0(h,rmask):
    """
    function rx0 = rx0(h,rmask)
 
    This function computes the bathymetry slope from a SCRUM NetCDF file.

    On Input:
       h           bathymetry at RHO-points.
       rmask       Land/Sea masking at RHO-points.
 
    On Output:
       rx0         Beckmann and Haidvogel grid stiffness ratios.
    """

    Lp, Mp = h.shape
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

    hx = np.zeros((L,Mp))
    hy = np.zeros((Lp,M))

    hx = abs(h[1:,:] - h[:-1,:]) / (h[1:,:] + h[:-1,:])
    hy = abs(h[:,1:] - h[:,:-1]) / (h[:,1:] + h[:,:-1])

    hx = hx * umask
    hy = hy * vmask

    rx0 = np.maximum(np.maximum(hx[:,:-1],hx[:,1:]), np.maximum(hy[:-1,:],hy[1:,:]))

    rmin = rx0.min()
    rmax = rx0.max()
    ravg = rx0.mean()
    rmed = np.median(rx0)

    print '  '
    print 'Minimum r-value = ', rmin
    print 'Maximum r-value = ', rmax
    print 'Mean    r-value = ', ravg
    print 'Median  r-value = ', rmed

    return rx0
