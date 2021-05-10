import RandSyst as RS
import System as S
from RandomParticleFunctions_v2 import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def Parabola(x,a,b,c):
    return a*x**2+b*x+c

uxxmin = -0.05
uxxmax = 0.05
NPoints = 100
SystemSize = 1
def GetEBulk(Mc,q0,check=False):
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
        plt.scatter(E[:,0],E[:,1])
    return E[np.argmin(E[:,1])]
def GetEBulk2(Mc,q0,check=False):
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
def GetL4MU(Mc=0, q0=0,check=False,Parameter = None):
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
            E = Sys.AffineDeformation(uxxmin,uxxmin)
        else :
            E = Sys.AffineDeformation(du, du)
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
def GetLambda(Mc=0,q0=0,check=False,Parameter = None):
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
            E = Sys.AffineDeformation(uxxmin,-uxxmin)
        else :
            E = Sys.AffineDeformation(du, -du)
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

def ComputePoissonRatio(Mc=0,q0=0,check=False,Parameter=None):
    L4MU = GetL4MU(Mc,q0,check,Parameter)
    Lambda = GetLambda(Mc,q0,check,Parameter)
    L = 0.5*(L4MU-Lambda)
    return L/(L4MU-L)
