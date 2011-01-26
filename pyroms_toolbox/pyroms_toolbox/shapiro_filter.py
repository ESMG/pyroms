import numpy as np

def shapiro1(Finp,order,scheme=1):
    '''
    This function applies a 1D shapiro filter to input 1D field.
    (ripped from Hernan G. Arango shapiro1.m)

     On Input:

         Finp        Field be filtered (1D array).
         order       Order of the Shapiro filter (2,4,8,16,...).
         scheme      Switch indicating the type of boundary scheme to use:
                         scheme = 1  =>  No change at wall, constant order.
                         scheme = 2  =>  Smoothing at wall, constant order.
                         scheme = 3  =>  No change at wall, reduced order.
                         scheme = 4  =>  Smoothing at wall, reduced order.
                         scheme = 5  =>  Periodic, constant order.
    '''

    fourk=np.array([2.500000e-1,   6.250000e-2,    1.562500e-2,    3.906250e-3, \
                    9.765625e-4,   2.44140625e-4,  6.103515625e-5, 1.5258789063e-5, \
                    3.814697e-6,   9.536743e-7,    2.384186e-7,    5.960464e-8, \
                    1.490116e-8,   3.725290e-9,    9.313226e-10,   2.328306e-10, \
                    5.820766e-11,  1.455192e-11,   3.637979e-12,   9.094947e-13])


    Im = Finp.shape[0]
    order2 = np.int(np.floor(order/2))

    cor=np.zeros((Im))
    Fcor=np.zeros((Im))

    #Compute filter correction.

    if (scheme == 1):
        #Scheme 1:  constant order and no change at wall.
        for n in range(order2):
            if n != order2:
                cor[0] = 2. * (Finp[0] - Finp[1])
                cor[Im-1] = 2. * (Finp[Im-1] - Finp[Im-2])
            else:
                cor[0] = 0.
                cor[Im-1] = 0.

            cor[1:-1] = 2. * Finp[1:-1] - Finp[0:-2] - Finp[2:]

        Fcor= cor * fourk[order2-1]

    elif scheme == 2:
        #Scheme 2:  constant order, smoothed at edges.
        for n in range(order2):
            cor[0] = 2. * (Finp[0] - Finp[1])
            cor[Im-1] = 2. * (Finp[Im-1] - Finp[Im-2])
            cor[1:-1] = 2. * Finp[1:-1] - Finp[0:-2] - Finp[2:]
            
        Fcor= cor * fourk[order2-1]

    elif scheme == 3:
        #Scheme 3:  reduced order and no change at wall.
        for n in range(order2):
            Istr = n
            Iend = Im-k+1
            if n == 1:
                cor[0] = 2. * (Finp[0] - Finp[1])
                cor[Im-1] = 2. * (Finp[Im-1] - Finp[Im-2])
                cor[1:-1] = 2. * Finp[1:-1] - Finp[0:-2] - Finp[2:]
            else:
                cor[Istr:Iend] = 2. * Finp[Istr:Iend] - Finp[Istr-1:Iend-1] - Finp[Istr+1:Iend+1]

            Fcor[Istr] = cor[Istr] * fourk[n]
            Fcor[Iend] = cor[Iend] * fourk[n]

        Fcor[0] = 0.
        Fcor[Istr:Iend] = cor[Istr:Iend] * fourk[order2]
        Fcor[Im] = 0.

    elif scheme == 4:
        #Scheme 4:  reduced order, smoothed at edges.
        for n in range(order2):
            Istr = n
            Iend = Im-k+1
            if n == 1:
                cor[0] = 2. * (Finp[0] - Finp[1])
                cor[Im-1] = 2. * (Finp[Im-1] - Finp[Im-2])
                cor[1:-1] = 2. * Finp[1:-1] - Finp[0:-2] - Finp[2:]
            else:
                cor[Istr:Iend] = 2. * Finp[Istr:Iend] - Finp[Istr-1:Iend-1] - Finp[Istr+1:Iend+1]

            Fcor[Istr] = cor[Istr] * fourk[n]
            Fcor[Iend] = cor[Iend] * fourk[n]

        Fcor[Istr:Iend] = cor[Istr:Iend] * fourk[order2]

    elif scheme == 5:
        #Scheme 5:  constant order, periodic.
        for n in range(order2):
            cor[0] = Finp[Im-2]
            cor[Im-1] = Finp[1]
            cor[1:-1] = 2. * Finp[1:-1] - Finp[0:-2] - Finp[2:]
            
        Fcor= cor * fourk[order2-1]

    #Apply correction.
    Fout = Finp - Fcor

    return Fout



def shapiro2(Finp,order,scheme=1,napp=1):
    '''
    This function applies a 2D shapiro filter to input 2D field.
    (ripped from Hernan G. Arango shapiro2.m)

     On Input:

         Finp        Field be filtered (2D array).
         order       Order of the Shapiro filter (2,4,8,16,...).
         scheme      Switch indicating the type of boundary scheme to use:
                         scheme = 1  =>  No change at wall, constant order.
                         scheme = 2  =>  Smoothing at wall, constant order.
                         scheme = 3  =>  No change at wall, reduced order.
                         scheme = 4  =>  Smoothing at wall, reduced order.
                         scheme = 5  =>  Periodic, constant order.
         napp        Number of Shapiro filter applications (optional). 
    '''

    Im, Jm = Finp.shape

    F=Finp.copy()
    Fout = np.zeros((Im, Jm))

    for n in range(napp):
  
        #Filter all rows.
        for j in range(Jm):
            Fraw = np.squeeze(F[:,j])
            Fraw = Fraw.T
            Fwrk = shapiro1(Fraw,order,scheme)
            Fout[:,j] = Fwrk.T

        #Filter all columns.
        for i in range(Im):
            Fraw = np.squeeze(Fout[i,:])
            Fwrk = shapiro1(Fraw,order,scheme)
            Fout[i,:] = Fwrk
  
    return Fout
