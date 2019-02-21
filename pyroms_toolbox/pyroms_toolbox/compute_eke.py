import numpy as np

def compute_eke(u, v):

    u = 100 * u #cm/s
    v = 100 * v #cm/s

    #move u and v to psi position
    u = 0.5 * (u[:,1:,:] + u[:,:-1,:])
    v = 0.5 * (v[:,:,1:] + v[:,:,:-1])

    u2 = u * u
    v2 = v * v

    uavg = u.mean(axis=0)
    vavg = v.mean(axis=0)

    u2avg = u2.mean(axis=0)
    v2avg = v2.mean(axis=0)

    eke = 0.5 * ((u2avg - uavg * uavg) + (v2avg - vavg * vavg))

    return eke

