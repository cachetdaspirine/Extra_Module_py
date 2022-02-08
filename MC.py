import numpy as np
import pathlib

CouplingRelation=np.load(str(pathlib.Path(__file__).parent.absolute())+'/CouplingRelation.npy',allow_pickle=True)
MainRelation=np.load(str(pathlib.Path(__file__).parent.absolute())+'/MainRelation.npy',allow_pickle=True)
VolumiqueRelationP=np.load(str(pathlib.Path(__file__).parent.absolute())+'/VolumiqueRelationP.npy',allow_pickle=True)
VolumiqueRelationM=np.load(str(pathlib.Path(__file__).parent.absolute())+'/VolumiqueRelationM.npy',allow_pickle=True)
#print(VolumiqueRelationP)
#print(VolumiqueRelationM)

def get_Mc(k=0,kc=0,eps=0,kA=0,kA2=False,Parameter = False):
    if Parameter:
        k,kc,eps,kA,kA2 = Parameter.k, Parameter.kc,Parameter.epsilon,Parameter.kA,Parameter.kA2
    if not kA2:
        kA2 = kA
    Mc = np.array([np.zeros(12,dtype=float) for _ in range(12)])
    Ap = (1+eps)**2*3**0.5/4.
    Am = (1-eps)**2*3**0.5/4.
    for R in CouplingRelation:
        Mc[R[0]]+=R[1]*kc/2.
    for R in MainRelation:
        Mc[R[0]]+=R[1]*k/2.
    for R in VolumiqueRelationP:
        Mc[R[0]]+=R[1]*kA/(2*Ap)
    for R in VolumiqueRelationM:
        Mc[R[0]]+=R[1]*kA2/(2*Am)

    Mc = (Mc+np.transpose(Mc))
    q0 = np.array([(1+eps)/(2*3**0.5),  #0
                    (1+eps)/2.,          #1
                    -(1-eps)/(2*3**0.5), #2
                    (1-eps)/2.,          #3
                    -(1+eps)/3**0.5,     #4
                    0.,                  #5
                    -(1-eps)/(2*3**0.5), #6
                    -(1-eps)/2.,         #7
                    (1+eps)/(2*3**0.5),  #8
                    -(1+eps)/2,          #9
                    (1-eps)/(3**0.5),    #10
                    0.])                 #11
    return Mc,q0

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
