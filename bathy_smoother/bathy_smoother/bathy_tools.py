import numpy as np

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def RoughnessMatrix(DEP, MSK):
    """
    RoughMat=GRID_RoughnessMatrix(DEP, MSK)

    ---DEP is the bathymetry of the grid
    ---MSK is the mask of the grid
    """

    eta_rho, xi_rho = DEP.shape

    Umat = np.array([[0, 1],
                    [1, 0],
                    [0, -1],
                    [-1, 0]])

    RoughMat = np.zeros(DEP.shape)

    for iEta in range(1,eta_rho-1):
        for iXi in range(1,xi_rho-1):
            if (MSK[iEta,iXi] == 1):
                rough = 0
                for i in range(4):
                    iEtaB = iEta + Umat[i,0]
                    iXiB = iXi + Umat[i,1]
                    if (MSK[iEtaB,iXiB] == 1):
                        dep1 = DEP[iEta,iXi]
                        dep2 = DEP[iEtaB,iXiB]
                        delta = abs((dep1 - dep2) / (dep1 + dep2))
                        rough = np.maximum(rough, delta)

                RoughMat[iEta,iXi] = rough


    return RoughMat
