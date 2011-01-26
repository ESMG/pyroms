def N2(rho, z, rho_0=1000.0):
    '''
    Return the stratification frequency
    
    Parameters
    ----------
    rho : array_like
        density [kg/m^3]
    z : array_like
        depths [m] (positive upward)
    
    Returns
    -------
    N2 : array_like
        Stratification frequency [1/s], where the vertical dimension (the
        first dimension) is one less than the input arrays.    
    '''
    rho = np.asarray(rho)
    z = np.asarray(z)
    assert rho.shape == z.shape, 'rho and z must have the same shape.'
    r_z = np.diff(rho, axis=0) / np.diff(z, axis=0)
    return -(9.8 / rho_0) * r_z

