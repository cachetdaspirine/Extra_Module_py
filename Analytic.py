import numpy as np
from scipy.special import iv

def Ff_Analytic(t,Gamma,nu):
    """
    return the free energy per volume unit for
    an infinit filament of width w.
    in the continuous unit:
    : Gamma : surface tension unit less
    : nu : Poisson ratio
    : t : width unit less
    : Return the energy in bulk free energy per volume unit
    """
    return Gamma*(1+nu)/t+1-(1+nu) * np.tanh(t/2)/t
def Ff_Analytic_Triangle(w,P):
    """
    Return the free energy per triangle for
    a set a of parameter:
    : w : width in number of triangle (2.w is the number of triangle)
    : P : an AnalyticToSimul or SimulToAnalytic object
    : Return the energy in bulk free energy per particle unit
    """
    if P.ParticleType == 'Triangle':
        t = np.sqrt(3)/2 * P.l0/P.l * w
    else :
        t = P.l0/P.l * w
    nu = P.nu
    Gamma = P.Gamma
    return P.FB*Ff_Analytic(t,Gamma,nu)
def FDisk_Anlytic(rhoSurf,rhoElastic,Gamma,nu):
    """
    return the free energy per volume unit for a Disk
    of radius rho in the continuous limit:
    : Gamma : surface tension unitless
    : nu : Poisson ratio
    : rhoSurf : radius unitless use for the surface free energy
    : rhoElastic : radius unitlesss use for the elastic free energy
    : Return the energy in bulk free energy per volume unit
    """
    return (Gamma*(1+nu)/rhoSurf+1-(1+nu)*((iv(0,rhoElastic)-iv(2,rhoElastic))/(iv(0,rhoElastic)\
                                    +iv(2,rhoElastic)+nu*(iv(0,rhoElastic)-iv(2,rhoElastic)))))
def FDisk_Analytic_Triangle(N,P,LowerBound = True):
    """
    return the free energy per triangle
    for a set of parameter:
    : N : is the number of triangle in the Hexagon
    : P : an AnalyticToSimul or SimulToAnalytic object
    : LowerBound : precise wether or not you want the LowerBound closest disk, or upperbound
    """
    if P.ParticleType == 'Triangle':
        rmin = np.sqrt(N)*P.l0/(P.l*2*np.sqrt(2))
        rmax = np.sqrt(N/6)*P.l0/P.l
        if LowerBound:
            return FDisk_Anlytic(rmin,rmin,P.Gamma,P.nu)*P.FB
        else :
            return FDisk_Anlytic(rmin,rmax,P.Gamma,P.nu)*P.FB
    else :
        r = 0.5+((9+12*(N-1))**0.5-3)/6.*P.l0/P.l
        return FDisk_Anlytic(r,r,P.Gamma,P.nu)*P.FB
