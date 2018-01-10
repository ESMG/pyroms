import numpy as np
from bathy_smoother import bathy_tools

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def smoothing_Positive_rx0(MSK, Hobs, rx0max):
    """
    This program use the direct iterative method from Martinho and Batteen (2006)
    The bathymetry is optimized for a given rx0 factor by increasing it.

    Usage:
    RetBathy = smoothing_Positive_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta,iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            LowerBound = RetBathy[iEtaN,iXiN] * (1-rx0max)/(1+rx0max)
                            if ((RetBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                RetBathy[iEta,iXi] = LowerBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    return RetBathy



def smoothing_Negative_rx0(MSK, Hobs, rx0max):
    """
    This program use an opposite methode to the direct iterative method from
    Martinho and Batteen (2006). This program optimizes the bathymetry for
    a given rx0 factor by decreasing it.

    Usage:
    RetBathy = smoothing_Negative_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            UpperBound = RetBathy[iEtaN, iXiN] * (1+rx0max)/(1-rx0max)
                            if (RetBathy[iEta,iXi] > (UpperBound + tol)):
                                IsFinished = 0
                                RetBathy[iEta, iXi] = UpperBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    return RetBathy



def smoothing_PositiveVolume_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This program use the direct iterative method from Martinho and Batteen (2006)
    The bathymetry is optimized for a given rx0 factor by increasing it. All depth
    are then multiplied by the coeficient K = Vol_init/Vol_final in order to
    insure volume conservation.

    Usage:
    RetBathy = smoothing_Positive_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho point
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    WorkBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            LowerBound = RetBathy[iEtaN, iXiN] * (1-rx0max)/(1+rx0max)
                            if ((WorkBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                WorkBathy[iEta, iXi] = LowerBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    VolOrig=0
    VolWork=0
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta, iXi] == 1):
                VolOrig = VolOrig + AreaMatrix[iEta,iXi] * Hobs[iEta,iXi]
                VolWork = VolWork + AreaMatrix[iEta,iXi] * WorkBathy[iEta,iXi]

    RetBathy = WorkBathy * (VolOrig / VolWork)

    return RetBathy



def smoothing_NegativeVolume_rx0(MSK, Hobs, rx0maxi, AreaMatrix):
    """
    This program use an opposite methode to the direct iterative method from
    Martinho and Batteen (2006). This program optimizes the bathymetry for
    a given rx0 factor by decreasing it. All depth are then multiplied by
    the coeficient K = Vol_init/Vol_final in order to insure volume conservation.

    Usage:
    RetBathy = smoothing_Negative_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho point
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    WorkBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            UpperBound = RetBathy[iEtaN, iXiN] * (1+rx0max)/(1-rx0max)
                            if (WorkBathy[iEta,iXi] > (UpperBound + tol)):
                                IsFinished = 0
                                WorkBathy[iEta, iXi] = UpperBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    VolOrig=0
    VolWork=0
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta, iXi] == 1):
                VolOrig = VolOrig + AreaMatrix[iEta,iXi] * Hobs[iEta,iXi]
                VolWork = VolWork + AreaMatrix[iEta,iXi] * WorkBathy[iEta,iXi]

    RetBathy = WorkBathy * (VolOrig / VolWork)

    return RetBathy



def smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This program use the Mellor-Ezer-Oey method (Mellor et al., 1994).
    The bathymetry is optimized for a given rx0 factor by doing a sequence
    of increase/decrease at adjacent cells.

    Usage:
    RetBathy, HmodifVal, ValueFct = smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho-points.
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    HmodifVal = 0
    TheMultiplier = (1 - rx0max) / (1 + rx0max)
    tol = 0.000001
    ValueFct = 0

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    Area = AreaMatrix[iEta, iXi]
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                            and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            AreaN = AreaMatrix[iEtaN,iXiN]
                            LowerBound = RetBathy[iEtaN,iXiN] * TheMultiplier
                            if ((RetBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                h = (TheMultiplier * RetBathy[iEtaN,iXiN] - RetBathy[iEta,iXi]) \
                                         / (AreaN + TheMultiplier * Area)
                                RetBathy[iEta,iXi] = RetBathy[iEta,iXi] + AreaN * h
                                RetBathy[iEtaN,iXiN] = RetBathy[iEtaN,iXiN] - Area * h
                                HmodifVal = HmodifVal + abs(h)
                                ValueFct = ValueFct + abs(h) * (Area + AreaN)

        if (IsFinished == 1):
            break

    H = AreaMatrix * Hobs * MSK
    TheBathymetry1 = H.sum()
    H = AreaMatrix * RetBathy * MSK
    TheBathymetry2 = H.sum()
    DeltaBathymetry = TheBathymetry1 - TheBathymetry2
    print('DeltaBathymetry = ', DeltaBathymetry)

    return RetBathy, HmodifVal, ValueFct


def smoothing_Laplacian_rx0(MSK, Hobs, rx0max):
    """
    This program use Laplacian filter.
    The bathymetry is optimized for a given rx0 factor by doing an iterated
    sequence of Laplacian filterings.

    Usage:
    RetBathy = smoothing_Laplacian_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    tol = 0.00001
    WeightMatrix = np.zeros((eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            WeightSum = 0
            for ineigh in range(4):
                iEtaN = iEta + ListNeigh[ineigh,0]
                iXiN = iXi + ListNeigh[ineigh,1]
                if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                      and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                    WeightSum = WeightSum + 1

            WeightMatrix[iEta,iXi] = WeightSum

    Iter = 1
    NumberDones = np.zeros((eta_rho, xi_rho))
    while(True):
        RoughMat = bathy_tools.RoughnessMatrix(RetBathy, MSK)
        Kbefore = np.where(RoughMat > rx0max)
        nbPtBefore = np.size(Kbefore, 1)
        realR = RoughMat.max()
        TheCorrect = np.zeros((eta_rho,xi_rho))
        IsFinished = 1
        nbPointMod = 0
        AdditionalDone = np.zeros((eta_rho, xi_rho))
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                Weight = 0
                WeightSum = 0
                for ineigh in range(4):
                    iEtaN = iEta + ListNeigh[ineigh,0]
                    iXiN = iXi + ListNeigh[ineigh,1]
                    if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                          and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                        Weight = Weight + RetBathy[iEtaN,iXiN]
                        AdditionalDone[iEtaN,iXiN] = AdditionalDone[iEtaN,iXiN] + NumberDones[iEta,iXi]

                TheWeight = WeightMatrix[iEta,iXi]
                WeDo = 0
                if TheWeight > tol:
                    if RoughMat[iEta,iXi] > rx0max:
                        WeDo = 1
                    if NumberDones[iEta,iXi] > 0:
                        WeDo = 1

                if WeDo == 1:
                    IsFinished = 0
                    TheDelta = (Weight - TheWeight * RetBathy[iEta,iXi]) / (2 * TheWeight)
                    TheCorrect[iEta,iXi] = TheCorrect[iEta,iXi] + TheDelta
                    nbPointMod = nbPointMod + 1
                    NumberDones[iEta,iXi] = 1

        NumberDones = NumberDones + AdditionalDone
        RetBathy = RetBathy + TheCorrect
        NewRoughMat = bathy_tools.RoughnessMatrix(RetBathy, MSK)
        Kafter = np.where(NewRoughMat > rx0max)
        nbPtAfter = np.size(Kafter, 1)
        TheProd = (RoughMat > rx0max) * (NewRoughMat > rx0max)
        nbPtInt = TheProd.sum()
        if (nbPtInt == nbPtAfter and nbPtBefore == nbPtAfter):
            eStr=' no erase'
        else:
            eStr='';
            NumberDones = np.zeros((eta_rho, xi_rho))

        print('Iteration #', Iter)
        print('current r=', realR, '  nbPointMod=', nbPointMod, eStr)
        print(' ')

        Iter = Iter + 1

        if (IsFinished == 1):
            break

    return RetBathy
