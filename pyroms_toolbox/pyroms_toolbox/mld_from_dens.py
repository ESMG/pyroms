import numpy as np
import pyroms

def mld_from_dens(dens,grd):

  z = grd.vgrid.z_r[0]

  surf_dens = dens[-1]
  mld_dens = surf_dens + 0.03 #threshold of 0.03 kg/m^3, following de Boyer Montegut et al., 2004

  mld, lon, lat = pyroms.tools.isoslice(z,dens,mld_dens, grd)

  return mld


