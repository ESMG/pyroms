import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from pyroms_toolbox import get_coast_line


def plot_coast_line(grd, proj=None, Cpos='rho'):
    '''
    plot_coast_line(grd, {proj})

    plot the coastline from the object grid. 
    proj=map (optional) is a Basemap object for
    projection.
    ''' 


    a = plt.gca()

    coast = get_coast_line(grd, Cpos=Cpos)
    c = np.array(coast)

    if proj is None:
        col = collections.LineCollection(c)
    else:
        cp = np.zeros(c.shape)
        for i in range(c.shape[0]):
            cp[i,:,0], cp[i,:,1] = proj(c[i,:,0], c[i,:,1])

        col = collections.LineCollection(cp)

 
    a.add_collection(col, autolim=True)
    col.set_color('k')
    #a.autoscale_view()
