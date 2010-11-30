# --- encoding: iso-8859-1 ---

"""
 Seawater -- Python functions for properties of sea water
    Bjørn Ådlandsvik <bjorn@imr.no>,
    Institute of Marine Research,
    Version 1.1, 13 November 2002

 Public functions:
 Density related
   dens(S,T,P)           Density of sea water              kg/m**3
   delta(S,T,P)          Specific volume anomaly           m**3/kg
   sigma(S,T,P)          Density anomaly                   kg/m**3
   drhodt(S,T,P)         Temperature derivative of density kg/(K*m**3)
   alpha(S,T,P)          Thermal expansion coefficient     1/K 
   drhods(S,T,P)         Salinity derivative of density    kg/m**3
   beta(S,T,P)           Salinity expansion coefficient
   
 Salinity related
   salt(R,T,P)           Salinity
   cond(S,T,P)           Conductivity ratio

 Heat related
   heatcap(S,T,P)        Heat capacity                     J/(kg*K)
   adtgrad(S,T,P)        Adiabatic lapse rate              K/dbar
   temppot(S,T,P,Pref)   Potential temperature             °C
   temppot0(S,T,P)       Potential temperature             °C
   
 Miscellaneous
   freezept(S,P)         Freezing point                    °C
   soundvel(S,T,P)       Sound velocity                    m/s
   depth(P,lat)          Depth                             m

 Arguments:
   S     = Salinity                 
   T     = Temperature               °C
   P     = Pressure                  dbar 
   R     = Conductivity ratio
   Pref  = Reference pressure        dbar
   lat   = Latitude                  deg

 References:
   [Bryden 1973], New polynomials for thermal expansion, adiabatic
   temperature gradient and potential temperature gradient of sea water
   Deep-Sea Res. 20, 401-408

   [UNESCO 1981], Tenth report of the joint panel on oceanographic
   tables and standards, Unesco technical papers in marine science 36.

   [UNESCO 1983], N.P. Fofonoff and R.C. Millard Jr., Algorithms for
   computation of fundamental properties of seawater, Unesco technical
   papers in marine science 44.

"""

# --- Exceptions ---
class OutOfRangeError(Exception): pass

from density import dens, svan, sigma, drhodt, alpha, drhods, beta
from salinity import salt, cond
from heat import heatcap, adtgrad, temppot, temppot0
from misc import freezept, soundvel, depth


