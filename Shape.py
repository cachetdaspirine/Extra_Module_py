import numpy as np
import copy

TopologieDownHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieUpHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieDownTriangle = [(1,0),(-1,0),(0,1)]
TopologieUpTriangle = [(1,0),(-1,0),(0,-1)]

def Parallel0(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<=size and j<i and j+i>=size/2 and i+j<=3*size/2-1:
                Res[i,j]=1
    return Res

def Parallel3(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<size and j<i and j+i>size/2 and i+j<=3*size/2-3:
                Res[i,j]=1
    return Res

def Parallel1(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<size and j<i and j+i>size/2-1 and i+j<=3*size/2-2:
                Res[i,j]=1
    return Res

def Parallel2(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<=size and j<i and j+i>size/2 and i+j<3*size/2-1:
                Res[i,j]=1
    return Res
def ParallelHex(size):
    Res=np.array([np.zeros(3*size,dtype=int) for _ in range(2*size)])
    Res[Res.shape[0]//2,Res.shape[1]//2] = 1
    for k in range(1,size):
        Resdt = copy.copy(Res)
        for i in range(Res.shape[0]):
            for j in range(Res.shape[1]):
                if Res[i,j]==1:
                    ij = (i,j)
                    for neigh in Get_Neighbors(Resdt,ij,Free = True,ParticleType = 'Hexagon'):
                        Resdt[neigh] = 1
        Res = Resdt

    return Res

#---------------------------------------------------------------
# I am not sure that it is necessary, but depending on the value
# of the size, there are 4 different way of  making  an  hexagon
# At the end the function Parallel ensure to return the  correct
# Hexagon
#---------------------------------------------------------------
#N = 6 * (Size/4)**2
def Star(Size,ParticleType = 'Triangle'):
    size = 4 * Size  # make sure size is pair
    Array = np.full((2*size-1,size +2 ),1)
    Array[:,0] = 0
    Array[:,-1] = 0
    for i in range(Array.shape[0]):
        for j in range(1,Array.shape[1]-1):
            if i < size//2  :
                if j-1 < size//2 :
                    if size//2-1-i > j-1:
                        Array[i,j] = 0
                else :
                    if j-1-size//2 > i :
                        Array[i,j] = 0
            elif i >= (size + size//2):
                if j-1 < size//2 :
                    if i - size - size//2 > j-1 -1 :
                        Array[i,j] = 0
                else :
                    if i-size - size//2 > size//2 -(j-1-size//2+2) :
                        Array[i,j] = 0
    Array[Array.shape[0]//2-size//2+3:Array.shape[0]//2+1+size//2-3,0] = 1
    Array[Array.shape[0]//2-size//2+3:Array.shape[0]//2+1+size//2-3,-1] = 1
    for i in range(Array.shape[0]):
        for j in range(1,Array.shape[1]-1):
            if i < size//2 -2 :
                if j-1 < size//2 :
                    if size//2-2-i == j-1:
                        Array[i,j] = 1
                    if size//2 -3-i == j-1 and i< size//2 -3:
                        Array[i,j] = 1
                else :
                    if j-2-size//2 <= i :
                        Array[i,j] = 1
                    if j-3-size//2 <= i and j<=size-1 :
                        Array[i,j] = 1
            elif i >= (size + size//2):
                if j-1 < size//2 :
                    if i - size - size//2 -1 <= j-1 -1 and i >size+size//2:
                        Array[i,j] = 1
                    if i - size - size//2 -2 <= j-1 -1 and i >size+size//2+1:
                        Array[i,j] = 1
                else :
                    if i-size//2-size < size+2-j and j < size :
                        Array[i,j] = 1
    return Array
def Parallel(size,ParticleType='Triangle'):
    if ParticleType=='Triangle':
        if size <4 :
            return np.array([[0,0,0],[0,1,0],[0,0,0]])
        if size%4==0:
            return Parallel0(size)
        elif size%4==1:
            return Parallel1(size)
        elif size%4==2:
            return Parallel2(size)
        elif size%4==3:
            return Parallel3(size)
    elif ParticleType=='Hexagon':
        return ParallelHex(size)
#---------------------------------------------------------------
# Given an array  of  0 and 1 : returns the number of 1, i.e the
# Number of particle
#---------------------------------------------------------------
def Np(Array):
    unique, counts = np.unique(Array, return_counts=True)
    return dict(zip(unique, counts))[1]

#---------------------------------------------------------------
# Make a fiber of certain width, and length.
#---------------------------------------------------------------
def Fiber(Width,Length,ParticleType='Triangle',type=1):
    if ParticleType=='Triangle':
        Res=np.array([np.zeros(Width+2,dtype=int) for _ in range(Length+2)])
        for i in range(1,Length+1):
            for j in range(1,Width+1):
                Res[i,j]=1
        return Res
    elif ParticleType=='Hexagon':
        if type == 1:
            return Fiber4(Width,Length)
        elif type == 2:
            return Fiber3(Width,Length)
        elif type == 3:
            return Fiber2(Width,Length)
def Fiber2(W,L):
    if W==1:
        return Fiber4(W,L)
    Res=np.array([np.zeros(L) for _ in range(W+int(L//2+0.5))])
    for l in range(L):
        if l%2==0:
            Res[l//2:W+l//2,l]=1
        else :
            Res[l//2+1:W+l//2,l]=1
    return np.flip(np.transpose(Res),0)
def Fiber3(W,L):
    if W==1:
        return Fiber4(W,L)
    Res=np.array([np.zeros(L) for _ in range(W+int(L//2+0.5))])
    for l in range(L):
        if l%2==0:
            Res[l//2:W+l//2,l]=1
        else :
            Res[l//2:W+l//2,l]=1
    return np.flip(np.transpose(Res),0)
def Fiber4(W,L):
    Res=np.array([np.zeros(W) for _ in range(L)])
    Res[:]=1
    return np.transpose(Res)
def D(ij,ij0,PType):
    if PType == 'Hexagon':
        return (((ij[0]-ij0[0])+0.5*((ij[1]-ij0[1])))**2+((ij[1]-ij0[1])*0.866)**2)**0.5
    else :
        return (((ij[0]-ij0[0])*0.5)**2+((ij[1]-ij0[1])*0.866)**2)**0.5
def Disk(L,ParticleType):
    if ParticleType == 'Hexagon':
        Ly,Lx = L,L
    else:
        Ly,Lx = L,L
    State = np.array([[1 for _ in range(Ly+1)] for _ in range(Lx+1)])
    ij0 = (Lx//2,Ly//2)
    R = Lx//4
    for i in range(State.shape[0]):
        for j in range(State.shape[1]):
            if D(ij0,(i,j),ParticleType) > R:
                State[i,j]=0
    if ParticleType == 'Triangle':
        for i in range(State.shape[0]):
            for j in range(State.shape[1]):
                if Get_Neighbors(State,(i,j),Occupied=True,ParticleType=ParticleType).__len__()== 0:
                    State[i,j] = 0
        AdjustDisk(State,ParticleType)
    return State
def AdjustDisk(State,ParticleType):
    for i in range(State.shape[0]):
        for j in range(State.shape[1]):
            if State[i,j]==0:
                if Get_Neighbors(State,(i,j),Occupied=True,ParticleType=ParticleType).__len__()>=2:
                    State[i,j]=1
def  Torrus(L,W,ParticleType):
    Ly,Lx = L,L
    State = np.full((Ly+1,Lx+1),1)
    #State = Parallel(L,ParticleType)
    #State = np.array([[1 for _ in range(Ly+1)] for _ in range(Lx+1)])
    ij0 = (Lx//2,Ly//2)
    RMax,RMin = Lx//4,Lx//4-W
    for i in range(State.shape[0]):
        for j in range(State.shape[1]):
            if D(ij0,(i,j),ParticleType) > RMax or D(ij0,(i,j),ParticleType) <= RMin:
            #if D(ij0,(i,j),ParticleType) <= RMin:
                State[i,j]=0
    return State
def Get_Neighbors(Array,ij,Occupied=False,Free=False,Border=False,ParticleType = 'Triangle'):
        Lx, Ly = Array.shape[0], Array.shape[1]
        if ParticleType == 'Triangle':
            TopologieUp = TopologieUpTriangle
            TopologieDown = TopologieDownTriangle
        elif ParticleType=='Hexagon':
            TopologieUp = TopologieUpHex
            TopologieDown = TopologieDownHex
        if (ij[0]+ij[1])%2==0:
            Res = np.array(TopologieDown)+np.array(ij)
        else :
            Res = np.array(TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        if not Border:
            Resreg=np.delete(Res,np.argwhere((Res[:,0]>=Lx) | (Res[:,0]<0) | (Res[:,1]>=Ly) | (Res[:,1]<0)),0)
        else:
            Resreg=Res
        #Build a numpy array of tuple
        Resbis=np.empty(Resreg.__len__(),dtype=object)
        Resbis[:] = list(zip(Resreg[:,0],Resreg[:,1]))
        if Border:
            Resbis = list(Resbis)
        #check the occupancie or not
        if Border :
            for n in reversed(range(Resbis.__len__())):
                if Resbis[n][0]>=0 and Resbis[n][1]>=0 and Resbis[n][0]<Lx and Resbis[n][1]<Ly:
                    if Array[Resbis[n]]!=0:
                        del Resbis[n]
        else :
            if Occupied:
                Resbis=Resbis[np.array([Array[r]==1 for r in Resbis])]
            elif Free:
                Resbis = Resbis[np.array([Array[r]==0 for r in Resbis])]

        return set(Resbis)
def SurfaceEnergy(Array,J=1.,ParticleType='Triangle'):
    surface = 0
    for i in range(Array.shape[0]):
        for j in range(Array.shape[1]):
            ij = (i,j)
            if Array[ij] == 1:
                surface+=Get_Neighbors(Array,ij,Free=True,Border=True,ParticleType=ParticleType).__len__()
    return J*surface
def GetSurfaceParticle(Array,ParticleType = 'Triangle'):
    SurfParticle = 0
    for i in range(Array.shape[0]):
        for j in range(Array.shape[1]):
            ij = (i,j)
            if Array[ij] == 1:
                if Get_Neighbors(Array,ij,Free=True,Border=True,ParticleType=ParticleType).__len__() != 0:
                    SurfParticle+=1
    return SurfParticle
def Lacunar(Size, Order =1, ParticleType='Triangle'):
    Array = Parallel(Size,ParticleType)
    if ParticleType == 'Hexagon':
        Ly,Lx = Array.shape[1]//2,Array.shape[0]//2
    else:
        Ly,Lx = Array.shape[1]//2,Array.shape[0]//2
    if Lx%2 == 1:
        Lx+=1
    if Ly%2 ==1 :
        Ly+=1
    # initialize the list of the particle to remove
    if ParticleType == 'Hexagon':
        To_Delete = {(Lx,Ly)}
    else :
        To_Delete = {(Lx-1,Ly),(Lx,Ly)}
    while To_Delete.__len__()!= 0:
        # Stor the initial To_delete
        #Old_To_Delete = To_Delete.copy()
        New_Delete = set()
        # add the neighboring ij's
        for ij in To_Delete:
            Array[ij] = 0
            for neigh in Get_Lacune(ij,Order,ParticleType):
                # if they follow two conditions :
                # They are in the box
                try :
                    # They are not already empty
                    if Array[neigh]==1:
                        New_Delete.add(neigh)
                except IndexError:
                    pass
        # Repete with the new added  one
        To_Delete = New_Delete
    return Array
def Get_Lacune(ij,Order,ParticleType):
    i,j=ij[0],ij[1]
    if Order ==0 and ParticleType == 'Hexagon':
        return  [(i-1,j+2),(i+1,j+1),(i+2,j-1),(i+1,j-2),(i-1,j-1),(i-2,j+1)]
    if Order == 1 and ParticleType == 'Hexagon':
        return [(j,j+2),(i+2,j),(i+2,j-2),(i,j-2),(i-2,j),(i-2,j+2)]
    if Order == 2 and ParticleType == 'Hexagon':
        return [(i-2,j+4),(i+2,j+2),(i+4,j-2),(i+2,j-4),(i-2,j-2),(i-4,j+2)]
    if ParticleType=='Triangle':
        return [(i+4,j),(i+2,j+2),(i-2,j+2),(i-4,j),(i-2,j-2),(i+2,j-2)]
