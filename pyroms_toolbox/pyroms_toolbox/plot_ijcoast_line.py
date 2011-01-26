import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from pyroms_toolbox import get_ijcoast_line


def plot_ijcoast_line(mask):
    '''
    plot_ijcoast_line(mask)

    plot the coastline from mask. 
    ''' 


    a = plt.gca()

    ijcoast = get_ijcoast_line(mask)
    c = np.array(ijcoast)

    col = collections.LineCollection(c)

    a.add_collection(col, autolim=True)
    col.set_color('k')
    #a.autoscale_view()
