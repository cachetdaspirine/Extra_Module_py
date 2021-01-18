import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
def delta(mu,nu):
    if mu==nu:
        return 1.
    else:
        return 0.
#This small routine take a list of position of center of mass of particle, compute the inertia matrix in 2D, then compute the principal axis and then the anisotropy of a list of particle. For a single aggregate!!
class Aggregate:
    def __init__(self,PositionFile=True,Array=True):
        if type(PositionFile)==str:
            Position=np.loadtxt(PositionFile,dtype=float,skiprows=1,delimiter=',')
        elif type(Array)==np.ndarray:
            Position=Array
        #In the normal file from the simulation X_G and Y_G are respectively the first and third column
        self.R=np.transpose(np.array([Position[:,0],Position[:,2]]))
        self.Npart=self.R.shape[0]
    def ComputeInertia(self):
        #start by computing the center of mass
        self.ComputeCenter()
        self.I=np.array([[0.,0.],[0.,0.]])
        for i in range(len(self.R)):
            for mu in range(2):
                for nu in range(2):
                    self.I[mu,nu]+=(-(self.R[i,mu]-self.R_G[mu])*(self.R[i,nu]-self.R_G[nu])+delta(mu,nu)*((self.R[i,0]-self.R_G[0])**2+(self.R[i,1]-self.R_G[1])**2))#/len(self.R)
    def ComputeCenter(self):
        self.R_G=np.zeros(2,dtype=float)
        self.R_G[0]=sum(self.R[:,0])/len(self.R)
        self.R_G[1]=sum(self.R[:,1])/len(self.R)
    def ComputePrincipale(self,ComputeInertiaMatrix=False):
        #Compute the inertia matrix if it doesn't exist or if we precise we want to compute it
        if not hasattr(self,'I') or not ComputeInertiaMatrix:
            self.ComputeInertia()
        self.w,self.v=LA.eigh(self.I)
        #The smallest eigen value means larger axis :
        if self.w[0]!=0 and self.w[1]!=0:
            self.Anisotropy=(1/self.w[0]-1/self.w[1])/(1/self.w[0]+1/self.w[1])
        else:
            self.Anisotropy=0
