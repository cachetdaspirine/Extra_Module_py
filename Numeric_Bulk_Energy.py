import numpy as np
import System as Sys
import RandSyst as RSys
import Shape as Sh
from MC import *

class BB:
    def __init__(self,OrderMax,Nmax,P,Expansion=False,Mc=False,q0=False):
        self.Nmax = Nmax
        self.Mc = Mc
        self.q0 = q0
        self.Expansion = Expansion
        self.OrderMax = OrderMax+1
        self.Systems = np.empty(self.OrderMax,dtype=object)
        if P.ParticleType=='Hexagon':
            for O in range(self.OrderMax):
                if O == 0 :
                    continue
                #self.Systems[O] = Sys.System(Sh.Lacunar(P.HSize(self.Nmax),O,P.ParticleType),Parameter = P)
                self.Systems[O] = self.MakeSystem(Sh.Lacunar(P.HSize(self.Nmax),O,P.ParticleType),P)
    def MakeSystem(self,State,P):
        if self.Expansion :
            if isinstance(self.Mc,np.ndarray) and isinstance(self.q0,np.ndarray):
                return RSys.System(self.Mc,self.q0,State)
            else :
                Mc,q0 = get_Mc(Parameter = P)
                return RSys.System(Mc,q0,State)
        else :
            return Sys.System(State = State,Parameter = P)
    def Get_E(self,Order,P):
        if Order == 0 and P.ParticleType=='Hexagon':
            return P.Flacune
        Array = Sh.Lacunar(P.HSize(self.Nmax),Order,P.ParticleType)
        Esurf = Sh.SurfaceEnergy(Array,P.J,P.ParticleType)
        Np = Sh.Np(Array)
        return (self.Systems[Order].Energy+Esurf)/Np

    def Get_Best_Bulk(self,P):
        if P.ParticleType == 'Triangle':
            return 0,np.inf
        E = [self.Get_E(O,P) for O in range(self.OrderMax)]
        return np.argmin(E),min(E)
