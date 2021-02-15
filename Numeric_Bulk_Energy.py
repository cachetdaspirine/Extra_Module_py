import numpy as np
import System as Sys
import Shape as Sh

class BB:
    def __init__(self,OrderMax,Nmax,P):
        self.Nmax = Nmax
        self.OrderMax = OrderMax+1
        self.Systems = np.empty(self.OrderMax,dtype=object)
        if P.ParticleType=='Hexagon':
            for O in range(self.OrderMax):
                if O == 0 :
                    continue
                self.Systems[O] = Sys.System(Sh.Lacunar(P.HSize(self.Nmax),O,P.ParticleType),Parameter = P)
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
