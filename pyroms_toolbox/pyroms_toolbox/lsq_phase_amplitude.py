import numpy as np
import matplotlib.pyplot as plt


def lsq_phase_amplitude(omega,ue,un,t):
    '''
    Amp, Pha = lsq_phase_amplitude(omega,ue,un,t)
    '''
    nc = omega.shape[0]
    m = 1 + 2*nc

    # Build Matrix
    c = np.ones((m,t.shape[0]))
    for i in range(nc):
        c[2*i+1,:]   =  np.cos(omega[i]*t)
        c[2*i+2,:] =  np.sin(omega[i]*t)

    A = np.zeros((m,m))
    for i in range(m):
        for j in range(m):
            A[i,j] = np.mean(c[i,:] * c[j,:])

    b1 = np.zeros(m)
    b2 = np.zeros(m)
    for i in range(m):
        b1[i] = np.mean(ue * c[i,:])
        b2[i] = np.mean(un * c[i,:])


    # Get solution 
    x1 = np.linalg.solve(A,b1)
    x2 = np.linalg.solve(A,b2)

    C = x1[1:2*nc:2]; D = x1[2:2*nc+1:2]

    Amp = np.sqrt(C*C + D*D)
    Pha = 180*np.arctan2(C,D)/np.pi

    for i in range(nc):
        if Pha[i] < 0:
            Pha[i] = Pha[i] + 360

    return Amp, Pha
