import RandSyst as RS
import System as S
from RandomParticleFunctions_v4 import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import copy
import Shape as Sh
import pathlib
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
    #V = Sys.Extension(0)*Sys.Extension(1)
    Sys.GetBulkEnergy()
    E0 = Sys.Energy
    for n in range(NPoints+1):
        uxx=(uxxmax-uxxmin)/NPoints*n+uxxmin
        if n ==0:
            E = Sys.AffineDeformation(uxxmin,uxxmin,Nodes)
        else :
            E = Sys.AffineDeformation(du, du,Nodes)
        l4mu.append([uxx,E])#(E-E0)/V
    l4mu = np.array(l4mu)
    p, conv = curve_fit(Parabola, l4mu[:, 0], l4mu[:, 1], p0=[0, 0, 0])
    if check:
        print(p)
        plt.plot(l4mu[:,0],l4mu[:,1],c='none')
        plt.scatter(l4mu[:,0],l4mu[:,1])
        plt.plot(l4mu[:,0],Parabola(l4mu[:,0],p[0],p[1],p[2]),label='Simulating particle l2mu')
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
        l.append([uxx,E])
    l = np.array(l)
    p, conv = curve_fit(Parabola, l[:, 0], l[:, 1], p0=[0, 0, 0])
    if check:
        print(p)
        plt.plot(l[:,0],l[:,1],c='none')
        plt.scatter(l[:,0],l[:,1])
        plt.plot(l[:,0],Parabola(l[:,0],p[0],p[1],p[2]),label='Simulating particle lambda')
    return p[0]

def ComputePoissonRatio(Mc=0,q0=0,check=False,Parameter=None,Nodes=[0,1]):
    L4MU = GetL4MU(Mc,q0,check,Parameter,Nodes)
    Lambda = GetLambda(Mc,q0,check,Parameter,Nodes)
    L = 0.5*(L4MU-Lambda)
    return L/(L4MU-L)
def compute_decoupled_nu(Mc,rho0,**kwargs):
    # This function measure the Poisson ratio of one of the two triangles
    # in the hexagonal particle. Taking a coupling matrix Mc, and a rest
    # configuration rho0 (rho0 corresponds to the 9 length) we construct
    # the according node position vector for a regular hexagon depending
    # on which triangle we wanna compute the Poisson ratio of. We Then
    # apply volumic, then pure shear deformation to measure nu. The def-
    # -formation is applied to the nodes position, and we compute the
    # length to get the energy accordingly.
    # Convert the length index from the Mert convention (see mapping of RandomParticleFunctions_v4)
    # To my mapping, see Thesis/

        #qregular_reduced = qregular
    #Measure dcoupled l2mu by applying a volumic deformation
    Mreduced,rho0reduced = reduce_coordinates(Mc,rho0,**kwargs)
    ux = np.array([1 for _ in range(6)])
    uy = np.array([1 for _ in range(6)])
    L2MU = deform_reduced_particle(Mreduced,qregular,rho0reduced,ux,uy,**kwargs)
    Lambda = deform_reduced_particle(Mreduced,qregular,rho0reduced,ux,-uy,**kwargs)
    L = 0.5*(L2MU-Lambda)
    return L/(L2MU-L)
def reduce_coordinates(Mc,rho0,**kwargs):
        IndexConversion = np.array([0,2,4,1,3,5,6,7,8])
        qregular = np.array([1/(2*3**.5),.5,-1/(2*3**.5),.5,-1/3**.5,0,-1/(2*3**0.5),-.5,1/(2*3**.5),-.5,1/(3**.5),0])
        if kwargs.get('triangle') == 'small':
            Mreduced = np.array([[Mc[j,i] for i in IndexConversion[3:6]] for j in IndexConversion[3:6]])
            rho0reduced = rho0[IndexConversion[3:6]]
            #qregular_reduced = qregular[[2,3,6,7,10,11]]
        elif kwargs.get('triangle')=='big':
            Mreduced = np.array([[Mc[i,j] for i in IndexConversion[0:3]] for j in IndexConversion[0:3]])
            #Mreduced = Mc[0:3,0:3]
            rho0reduced = rho0[IndexConversion[0:3]]
            #qregular_reduced = qregular[[0,1,4,5,8,9]]
        elif kwargs.get('triangle')=='both':
            Mreduced = np.array([[Mc[j,i] for i in IndexConversion[0:6]] for j in IndexConversion[0:6]])
            #Mreduced = Mc[0:6,0:6]
            rho0reduced = rho0[IndexConversion[0:6]]
            #qregular_reduced = qregular[[0,1,4,5,8,9,2,3,6,7,10,11]]
        else :#triangle == 'all':
            Mreduced = Mc
            rho0reduced = rho0
        return Mreduced,rho0reduced
def compute_stiffness_ratio(Mc,rho0,**kwargs):
    ux_fiber = np.load(str(pathlib.Path(__file__).parent.absolute())+'/ux_fiber.npy')
    uy_fiber = np.load(str(pathlib.Path(__file__).parent.absolute())+'/uy_fiber.npy')
    qregular = np.array([1/(2*3**.5),.5,-1/(2*3**.5),.5,-1/3**.5,0,-1/(2*3**0.5),-.5,1/(2*3**.5),-.5,1/(3**.5),0])
    if kwargs.get('eps'):
        eps=kwargs.get('eps')
    elif  kwargs.get('epsilon'):
        eps=kwargs.get('epsilon')
    else:
        eps=0.01
    q0 = get_Mc(eps=eps)[1]
    ux_bulk = qregular[0::2]-q0[0::2]
    uy_bulk = qregular[1::2]-q0[1::2]
    # if kwargs.get('triangle')=='small':
    #     ux_fiber = ux_fiber[[0,2,4]]
    #     uy_fiber = uy_fiber[[0,2,4]]
    #     ux_bulk = ux_bulk[[0,2,4]]
    #     uy_bulk = uy_bulk[[0,2,4]]
    # elif kwargs.get('triangle')=='big':
    #     ux_fiber = ux_fiber[[1,3,5]]
    #     uy_fiber = uy_fiber[[1,3,5]]
    #     ux_bulk = ux_bulk[[1,3,5]]
    #     uy_bulk = uy_bulk[[1,3,5]]
    ux_bulk,uy_bulk = ux_bulk/np.linalg.norm(ux_bulk),uy_bulk/np.linalg.norm(uy_bulk)
    Mreduced,rho0reduced = reduce_coordinates(Mc,rho0,**kwargs)
    return 0.5*deform_reduced_particle(Mreduced,qregular,rho0reduced,ux_bulk,uy_bulk,**kwargs)/deform_reduced_particle(Mreduced,qregular,rho0reduced,ux_fiber,uy_fiber,**kwargs)


def determine_local_coordinates(q,**kwargs):
    if kwargs.get('triangle')=='small':
        q = q[[2,3,6,7,10,11]]
        rho = np.array([(q[2*i]-q[(2*(i+1))%6])**2+(q[2*i+1]-q[(2*(i+1))%6+1])**2 for i in range(3)])
    elif kwargs.get('triangle')=='big':
        q = q[[0,1,4,5,8,9]]
        rho = np.array([(q[2*i]-q[(2*(i+1))%6])**2+(q[2*i+1]-q[(2*(i+1))%6+1])**2 for i in range(3)])
    elif kwargs.get('triangle')=='both':
        q = q[[0,1,4,5,8,9,2,3,6,7,10,11]]
        rho = np.array([(q[:6][2*i]-q[:6][(2*(i+1))%6])**2+(q[:6][2*i+1]-q[:6][(2*(i+1))%6+1])**2 for i in range(3)])
        rho = np.append(rho,np.array([(q[6:][2*i]-q[6:][(2*(i+1))%6])**2+(q[6:][2*i+1]-q[6:][(2*(i+1))%6+1])**2 for i in range(3)]))
    else:
        mapping = np.array([[0, 2, 4, 6, 8,10, 0, 4, 8],
                            [4, 6, 8,10, 0, 2, 2, 6,10]])
        rho = np.zeros(9)
        for ind_i in range(9):
            ind_1 = mapping[0, ind_i]
            ind_2 = mapping[1, ind_i]
            rho[ind_i] = (q[ind_1] - q[ind_2])**2+ (q[ind_1+1] - q[ind_2+1])**2
    return rho
def deform_reduced_particle(Mreduced,qregular_reduced,rho0reduced,ux,uy,**kwargs):
    def E(Mc,rho,rho0):
        return 0.5*np.dot((rho/rho0-rho0)/2.,np.dot(Mc,(rho/rho0-rho0)/2.))
    energy_list = np.zeros((NPoints,2))
    for n,amplitude in enumerate(np.linspace(uxxmin,uxxmax,NPoints)):
        # apply the deformation to the Nodes
        q_deformed = copy.copy(qregular_reduced)
        for i in range(ux.shape[0]):
            q_deformed[2*i] = qregular_reduced[2*i] * (1+ux[i]*amplitude)
            q_deformed[2*i+1] = qregular_reduced[2*i+1]*(1+uy[i]*amplitude)
        #compute the associated Length
        rho_deformed = determine_local_coordinates(q_deformed,**kwargs)
        energy_list[n] = [amplitude,E(Mreduced,rho_deformed,rho0reduced)]
    p, conv = curve_fit(Parabola, energy_list[:, 0], energy_list[:, 1], p0=[0, 0, 0])
    if kwargs.get('check'):
        print(p)
        plt.plot(energy_list[:,0],energy_list[:,1],c='none')
        plt.scatter(energy_list[:,0],energy_list[:,1])
        plt.plot(energy_list[:,0],Parabola(energy_list[:,0],p[0],p[1],p[2]),label='using matrix'+str(SignOfuxx*SignOfuyy))
    return p[0]
