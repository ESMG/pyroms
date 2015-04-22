import numpy as np

def get_littoral2(mask):
    '''
    littoral = get_littoral2(mask)

    return the littoral (first cells next to the shore) from 2D mask.
    Just like get_littoral, but without grabbing the diagonals.
    '''

    #get land point
    jidx, iidx = np.where(mask == 0)

    littoral = []

    for i in range(iidx.shape[0]):
        if jidx[i] != mask.shape[0]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]+1, iidx[i]]:
                littoral.append([jidx[i]+1, iidx[i]])

        if jidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]-1, iidx[i]]:
                littoral.append([jidx[i]-1, iidx[i]])

        if iidx[i] != mask.shape[1]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]+1]:
                littoral.append([jidx[i], iidx[i]+1])

        if iidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]-1]:
                littoral.append([jidx[i], iidx[i]-1])

    littoral = np.array(littoral)
    littoral = littoral.T[0], littoral.T[1]
    return littoral
