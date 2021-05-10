import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.sparse import csr_matrix
import sknetwork as skn
import math
def delta(mu,nu):
    if mu==nu:
        return 1.
    else:
        return 0.
def RDistance(IJ1,IJ2):
    return np.sqrt((np.mean(IJ1[0::2])-np.mean(IJ2[0::2]))**2+(np.mean(IJ1[1::2])-np.mean(IJ2[1::2]))**2)
def IsNeighbors(S1,S2):
    # S1 and S2 are two list of node position of size 12
    # we look at each nodes, if two are equal then they are neighbor
    # Convert the list of node into a set of tuple of (x,y)
    N1S = set(zip(truncate(S1[0::2],5),truncate(S1[1::2],5)))
    N2S = set(zip(truncate(S2[0::2],5),truncate(S2[1::2],5)))
    if len(N1S.intersection(N2S)) == 2:
        return True
    else:
        return False
@np.vectorize
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
#This small routine take a list of position of center of mass of particle, compute the inertia matrix in 2D, then compute the principal axis and then the anisotropy of a list of particle. For a single aggregate!!
class Aggregate:
    def __init__(self,PositionFile=True,Array=True):
        if type(PositionFile)==str:
            self.Position=np.loadtxt(PositionFile,dtype=float)
        elif type(Array)==np.ndarray:
            self.Position=Array
        self.MakeAdjacency()
        #In the normal file from the simulation X_G and Y_G are respectively the first and third column
        #self.R=np.transpose(np.array([Position[:,0],Position[:,2]]))
    def MakeAdjacency(self):
        Arr = np.zeros((self.Position.shape[0],self.Position.shape[0]),dtype=bool)
        for i in range(self.Position.shape[0]):
            for j in range(self.Position.shape[0]):
                if IsNeighbors(self.Position[i],self.Position[j]) or i==j:
                    Arr[i,j] = True
        self.A = csr_matrix(Arr)
    def ComputeShortestPaths(self):
        if not hasattr(self,'A'):
            self.MakeAdjacency()
        self.P = np.zeros((self.Position.shape[0],self.Position.shape[0]),dtype=object)
        self.LP = np.zeros((self.Position.shape[0],self.Position.shape[0]),dtype=int)
        for i1 in range(self.Position.shape[0]):
            for i2 in range(self.Position.shape[0]):
                self.P[i1,i2]=skn.path.shortest_path(self.A,i1,i2)
                self.LP[i1,i2] = len(self.P[i1,i2])
        #self.P = np.array([[skn.path.shortest_path(Adjacency,i1,i2) for i1 in range(self.Position.shape[0])] for i2 in self.Position.shape[0]//2])
        #self.P = self.P + self.P.T
    def TopologicalDistance(self,S1,S2):
        if not hasattr(self,'P'):
            self.ComputeShortestPaths()
        Distance = 0
        Path = self.P[S1,S2]
        for n in range(len(Path)-1):
            Distance+=RDistance(self.Position[Path[n]],self.Position[Path[n+1]])
        return Distance
    def XTopo(self,S):
        if not hasattr(self,'IG'):
            self.ComputeTopologicalCenter()
        XT = 0.
        Path = self.P[S,self.IG]
        for n in range(len(Path)-1):
            XT += abs(np.mean(self.Position[Path[n],0::2])-np.mean(self.Position[Path[n+1],0::2]))
        return XT
    def YTopo(self,S):
        if not hasattr(self,'IG'):
            self.ComputeTopologicalCenter()
        YT = 0.
        Path = self.P[S,self.IG]
        for n in range(len(Path)-1):
            YT += abs(np.mean(self.Position[Path[n],1::2])-np.mean(self.Position[Path[n+1],1::2]))
        return YT
    def ComputeTopologicalInertia(self):
        if not hasattr(self,'IG'):
            self.ComputeTopologicalCenter()
        self.IT=np.array([[0.,0.],[0.,0.]])
        for Site in range(self.Position.shape[0]):
            XY = (self.XTopo(Site),self.YTopo(Site))
            for mu in range(2):
                for nu in range(2):
                    self.IT[mu,nu] += - XY[mu]*XY[nu] + delta(mu,nu) * self.TopologicalDistance(self.IG,Site)**2
    def ComputeRealInertia(self):
        #start by computing the center of mass
        self.ComputeCenter()
        self.I=np.array([[0.,0.],[0.,0.]])
        for Site in self.Position:
            for mu in range(2):
                for nu in range(2):
                    self.I[mu,nu] += -(np.mean(Site[mu::2])-self.G[mu])*(np.mean(Site[nu::2])-self.G[nu])+delta(mu,nu)*((np.mean(Site[0::2])-self.G[0])**2 +(np.mean(Site[1::2])-self.G[1])**2)
                    #self.I[mu,nu]+=-(self.R[i,mu]-self.R_G[mu])*(self.R[i,nu]-self.R_G[nu])+delta(mu,nu)*((self.R[i,0]-self.R_G[0])**2+(self.R[i,1]-self.R_G[1])**2)#/len(self.R)
    def ComputeTopologicalCenter(self):
        if not hasattr(self,'P'):
            self.ComputeShortestPaths()
        self.IG = np.argmin(np.sum(self.LP,1))
    def ComputeCenter(self):
        self.G=np.zeros(2,dtype=float)
        self.G[0] = np.mean(self.Position[:,0::2])
        self.G[1] = np.mean(self.Position[:,1::2])
    def ComputePrincipale(self,Topology = False,ComputeInertiaMatrix=False):
        #Compute the inertia matrix if it doesn't exist or if we precise we want to compute it
        if Topology:
            if not hasattr(self,'IT') or not ComputeInertiaMatrix:
                self.ComputeTopologicalInertia()
            self.w,self.v=LA.eigh(self.IT)
        else:
            if not hasattr(self,'I') or not ComputeInertiaMatrix:
                self.ComputeRealInertia()
            self.w,self.v=LA.eigh(self.I)
        #The smallest eigen value means larger axis :
        if self.w[0]!=0 and self.w[1]!=0:
            self.Anisotropy=(1/self.w[0]-1/self.w[1])/(1/self.w[0]+1/self.w[1])
        else:
            self.Anisotropy=0
