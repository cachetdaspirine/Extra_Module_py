import RandSyst as RS
import System as S
from RandomParticleFunctions_v4 import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import copy
import Shape as Sh
from MC import *

def Parabola(x,a,b,c):
    return a*x**2+b*x+c
def line(x,a,b):
    return a*x+b

uxxmin = -0.05
uxxmax = 0.05
NPoints = 100
SystemSize = 1
def GetEBulk(Mc=False,q0=False,check=False,Parameter=False):
    if Parameter:
        Mc,q0 = get_Mc(Parameter=Parameter)
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

def ComputePoissonRatio(Mc=0,q0=0,check=False,Parameter=None,Nodes=[0,1]):
    L4MU = GetL4MU(Mc,q0,check,Parameter,Nodes)
    Lambda = GetLambda(Mc,q0,check,Parameter,Nodes)
    L = 0.5*(L4MU-Lambda)
    return L/(L4MU-L)
def compute_decoupled_nu(Mc,rho0,triangle='both'):
    # This function measure the Poisson ratio of one of the two triangles
    # in the hexagonal particle. Taking a coupling matrix Mc, and a rest
    # configuration rho0 (rho0 corresponds to the 9 length) we construct
    # the according node position vector for a regular hexagon depending
    # on which triangle we wanna compute the Poisson ratio of. We Then
    # apply volumic, then pure shear deformation to measure nu. The def-
    # -formation is applied to the nodes position, and we compute the
    # length to get the energy accordingly.
    qregular = np.array([(1)/(2*3**0.5),
                   (1)/2.,
                   -(1)/(2*3**0.5),
                   (1)/2.,
                   -(1)/3**0.5,
                   0.,
                   -(1)/(2*3**0.5),
                   -(1)/2.,
                   (1)/(2*3**0.5),
                   -(1)/2,
                   (1)/(3**0.5),
                   0.])
    if triangle == 'small':
        Mreduced = Mc[3:6,3:6]
        rho0reduced = rho0[3:6]
        qregular_reduced = qregular[[2,3,6,7,10,11]]
    elif triangle=='big':
        Mreduced = Mc[0:3,0:3]
        rho0reduced = rho0[0:3]
        qregular_reduced = qregular[[0,1,4,5,8,9]]
    elif triangle=='both':
        Mreduced = Mc[0:6,0:6]
        rho0reduced = rho0[0:6]
        qregular_reduced = qregular[[0,1,4,5,8,9,2,3,6,7,10,11]]
    #Measure dcoupled l2mu by applying a volumic deformation
    L2MU = deform_reduced_particle(Mreduced,qregular_reduced,rho0reduced,1,1)
    Lambda = deform_reduced_particle(Mreduced,qregular_reduced,rho0reduced,1,-1)
    L = 0.5*(L2MU-Lambda)
    return L/(L2MU-L)
def deform_reduced_particle(Mreduced,qregular_reduced,rho0reduced,SignOfuxx,SignOfuyy):
    def E(Mc,rho,rho0):
        return np.dot((rho-rho0),np.dot(Mc,rho-rho0))
    energy_list = np.zeros((NPoints,2))
    for n,uxx in enumerate(np.linspace(uxxmin,uxxmax,NPoints)):
        # apply the deformation to the Nodes
        q_deformed = copy.copy(qregular_reduced)
        q_deformed[0::2] = qregular_reduced[0::2] * (1+uxx*SignOfuxx)
        q_deformed[1::2] = qregular_reduced[1::2] * (1+uxx*SignOfuyy)
        #compute the associated Length
        if q_deformed.shape[0]//2 == 3:
            Nnodes = q_deformed.shape[0]//2
            rho_deformed = np.array([np.sqrt((q_deformed[2*(i%Nnodes)]-q_deformed[2*((i+1)%Nnodes)])**2+
                                            (q_deformed[2*(i%Nnodes)+1]-q_deformed[2*((i+1)%Nnodes)+1])**2)
                                            for i in range(Nnodes)])
        else:
            rho_deformed = np.zeros(6,dtype=float)
            Nnodes = q_deformed.shape[0]//4
            rho_deformed[:3] = np.array([np.sqrt((q_deformed[2*(i%Nnodes)]-q_deformed[2*((i+1)%Nnodes)])**2+
                                            (q_deformed[2*(i%Nnodes)+1]-q_deformed[2*((i+1)%Nnodes)+1])**2)
                                            for i in range(Nnodes)])
            rho_deformed[3:] = np.array([np.sqrt((q_deformed[2*(i%Nnodes+Nnodes)]-q_deformed[2*((i+1)%Nnodes+Nnodes)])**2+
                                            (q_deformed[2*(i%Nnodes+Nnodes)+1]-q_deformed[2*((i+1)%Nnodes+Nnodes)+1])**2)
                                            for i in range(Nnodes)])
        energy_list[n] = [uxx,E(Mreduced,rho_deformed,rho0reduced)]
    P, conv = curve_fit(Parabola, energy_list[:, 0], energy_list[:, 1], p0=[0, 0, 0])
    return P[0]
