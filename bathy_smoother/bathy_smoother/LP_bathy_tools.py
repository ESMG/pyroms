import numpy as np
from bathy_smoother import bathy_smoothing

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def GetIJS_rx0(MSK, DEP, r):

    eta_rho, xi_rho = DEP.shape
    print('eta_rho = ', eta_rho, '  xi_rho = ', xi_rho)

    nbVert = 0
    ListCoord = np.zeros((eta_rho, xi_rho))
    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1):
               nbVert = nbVert + 1
               ListCoord[iEta,iXi] = nbVert

    TotalNbVert = nbVert
    print('ListCoord built')
    print('Computing inequalities for r = ', r)

    TotalNbConstant = 0
    TotalNbEntry = 0
    for iEta in range(eta_rho-1):
       for iXi in range(xi_rho):
            if (MSK[iEta,iXi] == 1 and MSK[iEta+1,iXi] == 1):
                TotalNbConstant = TotalNbConstant + 2
                TotalNbEntry = TotalNbEntry + 4

    for iEta in range(eta_rho):
       for iXi in range(xi_rho-1):
            if (MSK[iEta,iXi] == 1 and MSK[iEta,iXi+1] == 1):
                TotalNbConstant = TotalNbConstant + 2
                TotalNbEntry = TotalNbEntry + 4

    TotalNbConstant = TotalNbConstant + 2 * TotalNbVert
    TotalNbEntry = TotalNbEntry + 4 * TotalNbVert

    Constant = np.zeros((TotalNbConstant,1))
    iList = np.zeros((TotalNbEntry,1))
    jList=np.zeros((TotalNbEntry,1))
    sList=np.zeros((TotalNbEntry,1))

    nbConst=0;
    nbEntry=0;
    for iEta in range(eta_rho-1):
       for iXi in range(xi_rho):
            if (MSK[iEta,iXi] == 1 and MSK[iEta+1,iXi] == 1):
                idx1 = ListCoord[iEta,iXi]
                idx2 = ListCoord[iEta+1,iXi]

                CST = (1+r) * DEP[iEta+1,iXi] + (-1+r) * DEP[iEta,iXi]
                Constant[nbConst,0] = CST
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx2
                sList[nbEntry,0] = -1-r
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx1
                sList[nbEntry,0] = 1-r
                nbEntry = nbEntry + 1

                CST = (1+r) * DEP[iEta,iXi] + (-1+r) * DEP[iEta+1,iXi]
                Constant[nbConst,0] = CST
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx1
                sList[nbEntry,0] = -r-1
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx2
                sList[nbEntry,0] = 1-r
                nbEntry = nbEntry + 1

    print('Inequalities for dh(iEta,iXi) and dh(iEta+1,iXi)')

    for iEta in range(eta_rho):
        for iXi in range(xi_rho-1):
            if (MSK[iEta,iXi] == 1 and MSK[iEta, iXi+1] == 1):
                idx1 = ListCoord[iEta,iXi]
                idx2 = ListCoord[iEta,iXi+1]

                CST = (1+r) * DEP[iEta,iXi+1] + (r-1) * DEP[iEta,iXi]
                Constant[nbConst,0] = CST
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx2
                sList[nbEntry,0] = -r-1
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx1
                sList[nbEntry,0] = 1-r
                nbEntry = nbEntry + 1

                CST = (1+r) * DEP[iEta,iXi] + (r-1) * DEP[iEta,iXi+1]
                Constant[nbConst,0] = CST
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx1
                sList[nbEntry,0] = -r-1
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx2
                sList[nbEntry,0] = 1-r
                nbEntry = nbEntry + 1

    print('Inequalities for dh(iEta,iXi) and dh(iEta,iXi+1)')

    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta,iXi] == 1):
                idx = ListCoord[iEta,iXi]

                Constant[nbConst,0] = 0
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = TotalNbVert + idx
                sList[nbEntry,0] = -1
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx
                sList[nbEntry,0] = 1
                nbEntry = nbEntry + 1

                Constant[nbConst,0] = 0
                nbConst = nbConst + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = TotalNbVert+idx
                sList[nbEntry,0] = -1
                nbEntry = nbEntry + 1
                iList[nbEntry,0] = nbConst
                jList[nbEntry,0] = idx
                sList[nbEntry,0] = -1
                nbEntry = nbEntry + 1

    print('Inequalities dh <= ad and -dh <= ad')

    print('rx0: nbEntry = ', nbEntry, '  nbConst = ', nbConst)
    print(' ')

    if (abs(nbEntry - TotalNbEntry) > 0):
        raise ValueError('We have a coding inconsistency for nbEntry. Please correct')

    if (abs(nbConst - TotalNbConstant) > 0):
        raise ValueError('We have a coding inconsistency for nbConst. Please correct')


    return iList, jList, sList, Constant



def GetIJS_maxamp(MSK, DEP, AmpConst):


    eta_rho, xi_rho = DEP.shape
    print('eta_rho = ', eta_rho, '  xi_rho = ', xi_rho)

    nbVert = 0
    ListCoord = np.zeros((eta_rho, xi_rho))
    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1):
               nbVert = nbVert + 1
               ListCoord[iEta,iXi] = nbVert

    TotalNbConstant = 0
    TotalNbEntry = 0
    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1):
                alpha = AmpConst[iEta,iXi]
                if (alpha < 9999):
                    TotalNbConstant = TotalNbConstant + 2
                    TotalNbEntry = TotalNbEntry + 2

    nbConst = 0
    nbEntry = 0
    Constant = np.zeros((TotalNbConstant,1))
    iList = np.zeros((TotalNbEntry,1))
    jList = np.zeros((TotalNbEntry,1))
    sList = np.zeros((TotalNbEntry,1))

    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1):
                idx = ListCoord[iEta,iXi]
                alpha = AmpConst[iEta,iXi]

                if (alpha < 9999):
                    Constant[nbConst,0] = alpha * DEP[iEta,iXi]
                    iList[nbEntry,0] = nbConst + 1
                    jList[nbEntry,0] = idx
                    sList[nbEntry,0] = -1
                    nbConst = nbConst + 1
                    nbEntry = nbEntry + 1

                    Constant[nbConst,0] = alpha * DEP[iEta,iXi]
                    iList[nbEntry,0] = nbConst + 1
                    jList[nbEntry,0] = idx
                    sList[nbEntry,0] = 1
                    nbConst = nbConst + 1
                    nbEntry = nbEntry + 1

    print('Inequalities |h^{new} - h^{old}| <= alpha h^{old}')
    print('maxamp: nbEntry = ', nbEntry, '  nbConst = ', nbConst)
    print(' ')

    if (abs(nbEntry - TotalNbEntry) > 0):
        raise ValueError('We have a coding inconsistency for nbEntry. Please correct')

    if (abs(nbConst - TotalNbConstant) > 0):
        raise ValueError('We have a coding inconsistency for nbConst. Please correct')


    return iList, jList, sList, Constant



def GetIJS_signs(MSK, SignConst):

    eta_rho, xi_rho = MSK.shape
    print('eta_rho = ', eta_rho, '  xi_rho = ', xi_rho)

    nbVert = 0
    ListCoord = np.zeros((eta_rho, xi_rho))
    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1):
               nbVert = nbVert + 1
               ListCoord[iEta,iXi] = nbVert

    TotalNbConstant = 0
    TotalNbEntry = 0
    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1 and SignConst[iEta,iXi] != 0):
                TotalNbConstant = TotalNbConstant + 1
                TotalNbEntry = TotalNbEntry + 1

    nbConst = 0
    nbEntry = 0
    Constant = np.zeros((TotalNbConstant,1))
    iList = np.zeros((TotalNbEntry,1))
    jList = np.zeros((TotalNbEntry,1))
    sList = np.zeros((TotalNbEntry,1))


    for iEta in range(eta_rho):
       for iXi in range(xi_rho):
           if (MSK[iEta,iXi] == 1 and SignConst[iEta,iXi] != 0):
               idx = ListCoord[iEta,iXi]

               Constant[nbConst,0] = 0
               nbConst = nbConst + 1
               iList[nbEntry,0] = nbConst
               jList[nbEntry,0] = idx
               if (SignConst[iEta,iXi] == 1):
                   sList[nbEntry,0] = -1
               elif (SignConst[iEta, iXi] == -1):
                   sList[nbEntry,0] = 1
               else:
                   raise ValueError('Wrong assigning please check SignConst')
               nbEntry = nbEntry + 1

    print('Inequalities dh >= 0 or dh <= 0')
    print('signs: nbEntry = ', nbEntry, '  nbConst = ', nbConst)
    print(' ')

    if (abs(nbEntry - TotalNbEntry) > 0):
        raise ValueError('We have a coding inconsistency for nbEntry. Please correct')

    if (abs(nbConst - TotalNbConstant) > 0):
        raise ValueError('We have a coding inconsistency for nbConst. Please correct')


    return iList, jList, sList, Constant



def MergeIJS_listings(iList1, jList1, sList1, Constant1, iList2, jList2, sList2, Constant2):

    # Suppose we have two sets of inequalities for two linear programs
    # with the same set of variables presented in sparse form.
    # The two descriptions are merge.

    nbConst1 = Constant1.shape[0]
    nbConst2 = Constant2.shape[0]
    nbEnt1 = iList1.shape[0]
    nbEnt2 = iList2.shape[0]

    Constant = np.zeros((nbConst1+nbConst2,1))
    iList = np.zeros((nbEnt1+nbEnt2,1))
    jList = np.zeros((nbEnt1+nbEnt2,1))
    sList = np.zeros((nbEnt1+nbEnt2,1))

    for iCons in range(nbConst1):
        Constant[iCons,0] = Constant1[iCons,0]

    for iCons in range(nbConst2):
        Constant[nbConst1+iCons,0] = Constant2[iCons,0]

    for iEnt in range(nbEnt1):
        iList[iEnt,0] = iList1[iEnt,0]
        jList[iEnt,0] = jList1[iEnt,0]
        sList[iEnt,0] = sList1[iEnt,0]

    for iEnt in range(nbEnt2):
        iList[nbEnt1+iEnt,0] = nbConst1 + iList2[iEnt,0]
        jList[nbEnt1+iEnt,0] = jList2[iEnt,0]
        sList[nbEnt1+iEnt,0] = sList2[iEnt,0]

    return iList, jList, sList, Constant


def GetBadPoints(MSK, DEP, rx0max):

    RetBathy = bathy_smoothing.smoothing_Positive_rx0(MSK, DEP, rx0max)
    K1 = np.where(RetBathy != DEP)

    eta_rho, xi_rho = MSK.shape
    MSKbad = np.zeros((eta_rho,xi_rho))
    MSKbad[K1] = 1

    return MSKbad


def Neighborhood(MSK, iEta, iXi, Kdist):

    eta_rho, xi_rho = MSK.shape
    MaxSiz = (2 * Kdist + 1) * (2 * Kdist + 1)
    ListNeigh = np.zeros((MaxSiz,2), dtype=np.int)
    ListStatus = -1 * np.ones((MaxSiz,1), dtype=np.int)
    ListKeys = np.zeros((MaxSiz,1), dtype=np.int)

    eKey = iEta + (eta_rho+1) * iXi
    ListNeigh[0,0] = iEta
    ListNeigh[0,1] = iXi
    ListStatus[0,0] = 0
    ListKeys[0,0] = eKey
    nbPt = 1

    List4dir = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    for iK in range(1,Kdist+1):
        nbPtOld = nbPt
        for iPt in range(nbPtOld):
            if (ListStatus[iPt,0] == iK-1):
                iEta = ListNeigh[iPt,0]
                iXi = ListNeigh[iPt,1]
                for ineigh in range(4):
                    iEtaN = iEta + List4dir[ineigh,0]
                    iXiN = iXi + List4dir[ineigh,1]
                    if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                            and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                        eKeyN = iEtaN + (eta_rho+1)*iXiN
                        Kf = np.where(ListKeys == eKeyN)
                        nbKf = np.size(Kf,1)
                        if (nbKf == 0):
                            ListNeigh[nbPt,0] = iEtaN
                            ListNeigh[nbPt,1] = iXiN
                            ListStatus[nbPt,0] = iK
                            ListKeys[nbPt,0] = eKeyN
                            nbPt = nbPt + 1

    ListNeighRet = ListNeigh[1:nbPt,:]

    return ListNeighRet


def ConnectedComponent(ListEdges, nbVert):
    """
    compute the vector of connected component belonging
    using a representation and an algorithm well suited
    for sparse graphs.
    """

    nbEdge = np.size(ListEdges, 0)
    ListDegree = np.zeros((nbVert,1), dtype=np.int)
    ListAdjacency = np.zeros((nbVert,10000), dtype=np.int)

    for iEdge in range(nbEdge):
        eVert = ListEdges[iEdge,0]
        fVert = ListEdges[iEdge,1]
        eDeg = ListDegree[eVert,0] + 1
        fDeg = ListDegree[fVert,0] + 1
        ListDegree[eVert,0] = eDeg
        ListDegree[fVert,0] = fDeg
        ListAdjacency[eVert,eDeg-1] = fVert
        ListAdjacency[fVert,fDeg-1] = eVert


    MaxDeg = ListDegree.max()
    ListAdjacency = ListAdjacency[:,:MaxDeg]

    ListVertexStatus = np.zeros((nbVert,1))
    ListHot = np.zeros((nbVert,1))
    ListNotDone = np.ones((nbVert,1))

    iComp = 0
    while(1):
        H = np.where(ListNotDone == 1)
        nb = np.size(H, 1)
        if (nb == 0):
            break;

        iComp = iComp + 1
        ListVertexStatus[H[0][0],0] = iComp
        ListHot[H[0][0],0] = 1
        while(1):
            H = np.where(ListHot == 1)
            ListNotDone[H] = 0
            ListNewHot = np.zeros((nbVert,1))
            for iH in range(np.size(H, 1)):
                eVert = H[0][iH]
                for iV in range(ListDegree[eVert, 0]):
                    ListNewHot[ListAdjacency[eVert, iV],0] = 1

            ListHot = ListNotDone * ListNewHot
            SumH = sum(ListHot)
            if (SumH == 0):
                break

            H2 = np.where(ListHot == 1)
            ListVertexStatus[H2] = iComp


    return ListVertexStatus
