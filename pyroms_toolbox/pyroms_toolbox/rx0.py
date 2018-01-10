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

    Mp, Lp = h.shape
    L=Lp-1
    M=Mp-1

    #  Land/Sea mask on U-points.
    umask = np.zeros((Mp,L))
    for j in range(Mp):
        for i in range(1,Lp):
            umask[j,i-1] = rmask[j,i] * rmask[j,i-1]

    #  Land/Sea mask on V-points.
    vmask = np.zeros((M,Lp))
    for j in range(1,Mp):
        for i in range(Lp):
            vmask[j-1,i] = rmask[j,i] * rmask[j-1,i]

    #-------------------------------------------------------------------
    #  Compute R-factor.
    #-------------------------------------------------------------------

    hx = np.zeros((Mp,L))
    hy = np.zeros((M,Lp))

    hx = abs(h[:,1:] - h[:,:-1]) / (h[:,1:] + h[:,:-1])
    hy = abs(h[1:,:] - h[:-1,:]) / (h[1:,:] + h[:-1,:])

    hx = hx * umask
    hy = hy * vmask

    rx0 = np.maximum(np.maximum(hx[:-1,:],hx[1:,:]),np.maximum(hy[:,:-1],hy[:,1:]))

    rmin = rx0.min()
    rmax = rx0.max()
    ravg = rx0.mean()
    rmed = np.median(rx0)

    print('  ')
    print('Minimum r-value = ', rmin)
    print('Maximum r-value = ', rmax)
    print('Mean    r-value = ', ravg)
    print('Median  r-value = ', rmed)

    return rx0
