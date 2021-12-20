import RandSyst as RS
import System as S
from RandomParticleFunctions_v4 import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import Shape as Sh

def Parabola(x,a,b,c):
    return a*x**2+b*x+c
def line(x,a,b):
    return a*x+b

uxxmin = -0.05
uxxmax = 0.05
NPoints = 100
SystemSize = 1
def GetEBulk(Mc,q0,check=False):
    #apply isotropic deformation
    State=np.full((SystemSize,SystemSize),1)
    Sys = RS.System(Mc, q0, State)
    du = (uxxmax-uxxmin)/NPoints
     #return Sys.GetBulkEnergy()
    E = [[0,Sys.GetBulkEnergy()]]
    for n in range(NPoints+1):
        uxx=(uxxmax-uxxmin)/NPoints*n+uxxmin
        if n ==0:
            E.append( [uxx,Sys.AffineDeformation(uxxmin,uxxmin)] )
        else :
            E.append( [uxx,Sys.AffineDeformation(du, du)] )
    E = np.array(E)
    if check:
        plt.plot(E[:,0],E[:,1],c='none')
        plt.scatter(E[:,0],E[:,1])
    #return uxx and the energy
    return E[np.argmin(E[:,1])]/Sh.Np(State)
def FindBestRegularHexagon(Mc,q0):
    qR =np.array([1/(2*3**0.5),(1)/2.,-(1)/(2*3**0.5),(1)/2.,-(1)/3**0.5,0,-(1)/(2*3**0.5),-(1)/2.,(1)/(2*3**0.5),-(1)/2,(1)/(3**0.5),0.])
    def E(Gamma):
        return 0.5*np.dot((1+Gamma)*qR-q0,np.dot(Mc,(1+Gamma)*qR-q0))
    res = minimize(E,0)
    return E(res.x),res.x
def GetEBulk2(Mc,q0,check=False):
    #apply shear deformation
    State=np.full((SystemSize,SystemSize),1)
    Sys = RS.System(Mc, q0, State)
    du = (uxxmax-uxxmin)/NPoints
     #return Sys.GetBulkEnergy()
    E = [[0,Sys.GetBulkEnergy()]]
    for n in range(NPoints+1):
        uxx=(uxxmax-uxxmin)/NPoints*n+uxxmin
        if n ==0:
            E.append( [uxx,Sys.AffineDeformation(uxxmin,-uxxmin)] )
        else :
            E.append( [uxx,Sys.AffineDeformation(du, -du)] )
    E = np.array(E)
    if check:
        plt.scatter(E[:,0],E[:,1])
    return E[np.argmin(E[:,1])]
def GetRealBulk(Mc,q0,*arg,**kwargs):
    Size = [10,15,20,30,40]
    E = list()
    for s in Size:
        NP = Sh.Np(Sh.Parallel(s,ParticleType='Hexagon'))
        E.append([NP,RS.System(Mc,q0,Sh.Parallel(s,ParticleType='Hexagon')).Energy/NP])
    E = np.array(E)
    p, conv = curve_fit(line, 1/np.sqrt(E[:, 0]), E[:, 1], p0=[0, 0])
    if 'Check' in kwargs or 'check' in kwargs:
        if kwargs.get('Check') or kwargs.get('check'):
            plt.plot(1/E[:,0]**0.5,line(1/E[:,0]**0.5,p[0],p[1]))
            plt.scatter(1/E[:,0]**0.5,E[:,1])
    return p[1]


def GetL4MU(Mc=0, q0=0,check=False,Parameter = None,Nodes=[0,1]):
    # Computation variables
    l4mu = list()
    du = (uxxmax-uxxmin)/NPoints
    State=np.full((SystemSize,SystemSize),1)
    ##########################################
    if Parameter :
        Sys = S.System(State,Parameter=Parameter)
    else:
        Sys = RS.System(Mc, q0, State)
    V = Sys.Extension(0)*Sys.Extension(1)
    Sys.GetBulkEnergy()
    E0 = Sys.Energy
    for n in range(NPoints+1):
        uxx=(uxxmax-uxxmin)/NPoints*n+uxxmin
        if n ==0:
            E = Sys.AffineDeformation(uxxmin,uxxmin,Nodes)
        else :
            E = Sys.AffineDeformation(du, du,Nodes)
        l4mu.append([uxx,(E-E0) / V])
    l4mu = np.array(l4mu)
    p, conv = curve_fit(Parabola, l4mu[:, 0], l4mu[:, 1], p0=[0, 0, 0])
    if check:
        fig = plt.figure(figsize=(8,6))
        print(p)
        plt.plot(l4mu[:,0],l4mu[:,1],c='none')
        plt.scatter(l4mu[:,0],l4mu[:,1])
        plt.plot(l4mu[:,0],Parabola(l4mu[:,0],p[0],p[1],p[2]))
    return p[0]
def GetLambda(Mc=0,q0=0,check=False,Parameter = None,Nodes=[0,1]):
    # Computation variables
    l = list()
    du = (uxxmax-uxxmin)/NPoints
    State=np.full((SystemSize,SystemSize),1)
    ##########################################
    if Parameter :
        Sys = S.System(State,Parameter=Parameter)
    else:
        Sys = RS.System(Mc, q0, State)
    V = Sys.Extension(0)*Sys.Extension(1)
    Sys.GetBulkEnergy()
    E0 = Sys.Energy
    for n in range(NPoints+1):
        uxx=(uxxmax-uxxmin)/NPoints*n+uxxmin
        if n ==0:
            E = Sys.AffineDeformation(uxxmin,-uxxmin,Nodes)
        else :
            E = Sys.AffineDeformation(du, -du,Nodes)
        l.append([uxx,(E-E0) / V])
    l = np.array(l)
    p, conv = curve_fit(Parabola, l[:, 0], l[:, 1], p0=[0, 0, 0])
    if check:
        fig = plt.figure(figsize=(8,6))
        print(p)
        plt.plot(l[:,0],l[:,1],c='none')
        plt.scatter(l[:,0],l[:,1])
        plt.plot(l[:,0],Parabola(l[:,0],p[0],p[1],p[2]))
    return p[0]

def ComputePoissonRatio(Mc=0,q0=0,check=False,Parameter=None,Nodes = [0,1]):
    L4MU = GetL4MU(Mc,q0,check,Parameter,Nodes)
    Lambda = GetLambda(Mc,q0,check,Parameter,Nodes)
    L = 0.5*(L4MU-Lambda)
    return L/(L4MU-L)
