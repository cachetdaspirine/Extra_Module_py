import math
from scipy.optimize import minimize
from scipy.optimize import brentq,newton
import numpy as np
import MeasurePoisson as MP
import random as rd
import matplotlib.pyplot as plt
from RandomParticleFunctions_v4 import *

class Mcq0:
    def __init__(self,seed=None,length=None,distribution='uniform',*argv,**kwargs):
        if seed:
            self.seed = seed
            np.random.seed(seed)
            rd.seed(seed)
        else :
            self.seed = np.random.randint(0,100000000)
        if length :
            self.Mc, self.rho0, self.e1,self.e2 = RandomParticle(self.seed,length,distribution=distribution)
            self.length=length
            self.vals,self.vect = np.linalg.eigh(self.Mc)
        else :
            self.length = 0.
            self.Mc, self.rho0, self.e1,self.e2 = RandomMatrix(self.seed,distribution=distribution)
            self.vals,self.vect = np.linalg.eigh(self.Mc)
        if kwargs.get('CheckEigen'):
            eigen = list()
            for p in np.arange(-1,1.,0.1):
                Mc,q0,e1,e2 = RandomParticle(self.R,p)
                arr = np.append(np.asarray(FindEigenValues(Mc)),p)
                if arr[0]>-10**-10:
                    eigen.append(arr[3])
                else :
                    eigen.append(arr[0])
            eigen = np.array(eigen)
            plt.plot(np.arange(-1,1.,0.1),eigen)
        if kwargs.get('GetPoisson'):
            self.FB,Gamma = MP.FindBestRegularHexagon(self.Mc,self.q0)
            self.nu = MP.ComputePoissonRatio(self.Mc,self.q0)
        if kwargs.get('CheckBulk'):
            self.RealFB = MP.GetRealBulk(self.Mc,self.q0)
    def GetSteadyParticle(self,seed=None,pressure=0.,percent=1):
        if seed:
            rd.seed(seed)
        self.R = rd.randint(0,1000000)
        Mc,q0,e1,e2 = RandomParticle(self.R,pressure,percent=percent)
        while FindEigenValues(Mc)[0]<-10**(-10):
            self.R = rd.randint(0,1000000)
            Mc,q0,e1,e2 = RandomParticle(self.R,pressure,percent=percent)
        return Mc,q0,e1,e2
