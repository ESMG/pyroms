import numpy as np

def get_ijcoast_line(mask):
    '''
    ijcoast = get_ijcoast_line(mask)

    return the ij coastline from mask 
    '''

    #get land point
    jidx, iidx = np.where(mask == 0)

    ijcoast = []

    for i in range(iidx.shape[0]):
        if jidx[i] != mask.shape[0]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]+1, iidx[i]]:
                ic = ([iidx[i]-0.5, iidx[i]+1-0.5])
                jc = ([jidx[i]+1-0.5, jidx[i]+1-0.5])
                seg = list(zip(ic,jc))
                ijcoast.append(seg)

        if jidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i]-1, iidx[i]]:
                ic = ([iidx[i]-0.5, iidx[i]+1-0.5])
                jc = ([jidx[i]-0.5, jidx[i]-0.5])
                seg = list(zip(ic,jc))
                ijcoast.append(seg)

        if iidx[i] != mask.shape[1]-1:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]+1]:
                ic = ([iidx[i]+1-0.5, iidx[i]+1-0.5])
                jc = ([jidx[i]-0.5, jidx[i]+1-0.5])
                seg = list(zip(ic,jc))
                ijcoast.append(seg)

        if iidx[i] != 0:
            if mask[jidx[i], iidx[i]] != mask[jidx[i], iidx[i]-1]:
                ic = ([iidx[i]-0.5, iidx[i]-0.5])
                jc = ([jidx[i]-0.5, jidx[i]+1-0.5])
                seg = list(zip(ic,jc))
                ijcoast.append(seg)

    return np.array(ijcoast)
