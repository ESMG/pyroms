import numpy as np
from bathy_smoother import LP_bathy_tools
from bathy_smoother import LP_tools
from bathy_smoother import bathy_tools

import matplotlib.pyplot as plt

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def LP_smoothing_rx0(MSK, Hobs, rx0max, SignConst, AmpConst):
    """
    This program perform a linear programming method in order to 
    optimize the bathymetry for a fixed factor r.
    The inequality |H(e)-H(e')| / (H(e)-H(e')) <= r where H(e)=h(e)+dh(e) 
    can be rewritten as two linear inequalities on dh(e) and dh(e'). 
    The optimal bathymetry is obtain by minimising the perturbation
    P = sum_e(|dh(e)| under the above inequalitie constraintes.

    Usage:
    NewBathy = LP_smoothing_rx0(MSK, Hobs, rx0max, SignConst, AmpConst)
   
    ---MSK(eta_rho,xi_rho) is the mask of the grd
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---SignConst(eta_rho,xi_rho) matrix of 0, +1, -1
         +1  only bathymetry increase are allowed.
         -1  only bathymetry decrease are allowed.
         0   increase and decrease are allowed.
         (put 0 if you are indifferent)
    ---AmpConst(eta_rho,xi_rho)  matrix of reals.
         coefficient alpha such that the new bathymetry should
         satisfy to  |h^{new} - h^{raw}| <= alpha h^{raw}
         (put 10000 if you are indifferent)
    """

    eta_rho, xi_rho = MSK.shape

    iList, jList, sList, Constant = LP_bathy_tools.GetIJS_rx0(MSK, Hobs, rx0max)

    iListApp, jListApp, sListApp, ConstantApp = LP_bathy_tools.GetIJS_maxamp(MSK, Hobs, AmpConst)

    iList, jList, sList, Constant = LP_bathy_tools.MergeIJS_listings(iList, jList, sList, Constant, iListApp, jListApp, sListApp, ConstantApp)

    iListApp, jListApp, sListApp, ConstantApp = LP_bathy_tools.GetIJS_signs(MSK, SignConst)

    iList, jList, sList, Constant = LP_bathy_tools.MergeIJS_listings(iList, jList, sList, Constant, iListApp, jListApp, sListApp, ConstantApp)

    TotalNbVert = int(MSK.sum())

    ObjectiveFct = np.zeros((2*TotalNbVert,1))
    for iVert in range(TotalNbVert):
        ObjectiveFct[TotalNbVert+iVert,0] = 1

    ValueFct, ValueVar, testfeasibility = LP_tools.SolveLinearProgram(iList, jList, sList, Constant, ObjectiveFct)
    if (testfeasibility == 0):
        NewBathy = NaN * np.ones((eta_rho,xi_rho))
        raise ValueError('Feasibility test failed. testfeasibility = 0.')

    correctionBathy = np.zeros((eta_rho,xi_rho))
    nbVert = 0
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta,iXi] == 1):
                correctionBathy[iEta,iXi] = ValueVar[nbVert]
                nbVert = nbVert + 1

    NewBathy = Hobs + correctionBathy
    RMat = bathy_tools.RoughnessMatrix(NewBathy, MSK)
    MaxRx0 = RMat.max()
    print('rx0max = ', rx0max, '  MaxRx0 = ', MaxRx0)

    return NewBathy





def LP_smoothing_rx0_heuristic(MSK, Hobs, rx0max, SignConst, AmpConst):
    """
    This program perform a linear programming method in order to 
    optimize the bathymetry for a fixed factor r.
    The inequality |H(e)-H(e')| / (H(e)-H(e')) <= r where H(e)=h(e)+dh(e) 
    can be rewritten as two linear inequalities on dh(e) and dh(e'). 
    The optimal bathymetry is obtain by minimising the perturbation
    P = sum_e(|dh(e)| under the above inequalitie constraintes.
    In order to reduce the computation time, an heurastic method is
    used.

    Usage:
    NewBathy = LP_smoothing_rx0_heuristic(MSK, Hobs, rx0max, SignConst, AmpConst)
    
    ---MSK(eta_rho,xi_rho) is the mask of the grd
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---SignConst(eta_rho,xi_rho) matrix of 0, +1, -1
         +1  only bathymetry increase are allowed.
         -1  only bathymetry decrease are allowed.
         0   increase and decrease are allowed.
         (put 0 if you are indifferent)
    ---AmpConst(eta_rho,xi_rho)  matrix of reals.
         coefficient alpha such that the new bathymetry should
         satisfy to  |h^{new} - h^{raw}| <= alpha h^{raw}
         (put 10000 if you are indifferent)
    """


    # the points that need to be modified
    MSKbad = LP_bathy_tools.GetBadPoints(MSK, Hobs, rx0max)

    eta_rho, xi_rho = MSK.shape

    Kdist = 5

    Kbad = np.where(MSKbad == 1)
    nbKbad = np.size(Kbad,1)
    ListIdx = np.zeros((eta_rho,xi_rho), dtype=np.int)
    ListIdx[Kbad] = list(range(nbKbad))

    ListEdges = []
    nbEdge = 0
    for iK in range(nbKbad):
        iEta, iXi = Kbad[0][iK], Kbad[1][iK]
        ListNeigh = LP_bathy_tools.Neighborhood(MSK, iEta, iXi, 2*Kdist+1)
        nbNeigh = np.size(ListNeigh, 0)
        for iNeigh in range(nbNeigh):
            iEtaN, iXiN = ListNeigh[iNeigh]
            if (MSKbad[iEtaN,iXiN] == 1):
                idx = ListIdx[iEtaN,iXiN]
                if (idx > iK):
                    nbEdge = nbEdge + 1
                    ListEdges.append([iK, idx])

    ListEdges = np.array(ListEdges)
    ListVertexStatus = LP_bathy_tools.ConnectedComponent(ListEdges, nbKbad)
    nbColor = ListVertexStatus.max()

    NewBathy = Hobs.copy()
    for iColor in range(1,nbColor+1):
        print('---------------------------------------------------------------')
        MSKcolor = np.zeros((eta_rho, xi_rho))
        K = np.where(ListVertexStatus == iColor)
        nbK = np.size(K,1)
        print('iColor = ', iColor, '  nbK = ', nbK)
        for iVertex in range(nbKbad):
            if (ListVertexStatus[iVertex,0] == iColor):
                iEta, iXi = Kbad[0][iVertex], Kbad[1][iVertex]
                MSKcolor[iEta, iXi] = 1
                ListNeigh = LP_bathy_tools.Neighborhood(MSK, iEta, iXi, Kdist)
                nbNeigh = np.size(ListNeigh, 0)
                for iNeigh in range(nbNeigh):
                    iEtaN, iXiN = ListNeigh[iNeigh]
                    MSKcolor[iEtaN,iXiN] = 1
        K = np.where(MSKcolor == 1)
        MSKHobs = np.zeros((eta_rho, xi_rho))
        MSKHobs[K] = Hobs[K].copy()
        TheNewBathy = LP_smoothing_rx0(MSKcolor, MSKHobs, rx0max, SignConst, AmpConst)
        NewBathy[K] = TheNewBathy[K].copy()

    print('Final obtained bathymetry')
    RMat = bathy_tools.RoughnessMatrix(NewBathy, MSK)
    MaxRx0 = RMat.max()
    print('rx0max = ', rx0max, '  MaxRx0 = ', MaxRx0)

    return NewBathy

