import numpy as np
import System as Sys
import RandSyst as RSys
import Shape as Sh
from MC import *


class BD:
    def __init__(self,Nmax,P,Expansion = False,Mc=False,q0=False):
        self.Expansion = Expansion
        self.Mc = Mc
        self.q0 = q0
        self.Nmax = Nmax
        self.Range = P.HRange(Nmax)
        self.NList = np.array([Sh.Np(Sh.Parallel(r,P.ParticleType)) for r in self.Range])
        #self.Systems = np.empty(self.Range.shape[0],dtype=object)
        self.SystemEnergy = np.zeros(self.Range.shape[0],dtype=float)
        #self.Systems = np.array([
        #Sys.System(
        #Sh.Parallel(r,ParticleType=P.ParticleType),
        #eps=P.epsilon,
        #Kcoupling=P.kc,
        #Kvol=P.kA,
        #ParticleType=P.ParticleType)
        #for r in self.Range])
    def MakeSystem(self,State,P):
        if self.Expansion :
            if self.Mc and self.q0:
                return RSys.System(self.Mc,self.q0,State)
            else:
                Mc,q0 = get_Mc(Parameter = P)
                return RSys.System(Mc,q0,State)
        else :
            return Sys.System(State = State,Parameter = P)
    def Get_E(self,n,P):
        if n>= self.Range.shape[0]:
            n=-1
        DiskArray = Sh.Parallel(self.Range[n],ParticleType=P.ParticleType)
        #if self.Systems[n]:
        #    return (self.Systems[n].Energy+Sh.SurfaceEnergy(DiskArray,J=P.J,ParticleType=P.ParticleType))/Sh.Np(DiskArray)
        if self.SystemEnergy[n]:
            return (self.SystemEnergy[n]+Sh.SurfaceEnergy(DiskArray,J=P.J,ParticleType=P.ParticleType))/Sh.Np(DiskArray)
        else :
            #self.SystemEnergy[n] =Sys.System(
            #DiskArray,
            #eps=P.epsilon,
            #Kcoupling=P.kc,
            #Kvol=P.kA,
            #ParticleType=P.ParticleType).Energy
            self.SystemEnergy[n] = self.MakeSystem(DiskArray,P).Energy
            return (self.SystemEnergy[n]+Sh.SurfaceEnergy(DiskArray,J=P.J,ParticleType=P.ParticleType))/Sh.Np(DiskArray)
    def Get_Best_Disk(self,P):
        n1 = 0
        E1 = self.Get_E(n1,P)
        n2 = 1
        E2 = self.Get_E(n2,P)
        while E2<E1:
            n1=n2
            n2+=1
            E1=E2
            E2 = self.Get_E(n2,P)
            if n2==self.Range.shape[0]-1:
                return Sh.Np(Sh.Parallel(self.Range[n2],ParticleType=P.ParticleType)),E2
        return Sh.Np(Sh.Parallel(self.Range[n1],ParticleType=P.ParticleType)),E1
