import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.sparse import csr_matrix
import sknetwork as skn
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict
import random
cdict = {'blue':   ((0.0,  0.9, 0.9),
                    (0.5,  0.4, 0.4),
                    (1.0,  0.1, 0.1)),

         'green': ((0.0,  0.5, 0.5),
                   (0.5, 1, 1),
                   (1.0,  0.3, 0.3)),

         'alpha': ((0.0,  1, 1),
                   (0.5, 0.8, 0.8),
                   (1.0,  1, 1)),

         'red':  ((0.0,  0.4, 0.4),
                  (0.5,  0.5, 0.5),
                  (1.0,  0.9, 0.9)),
         }
cm = LinearSegmentedColormap('my_colormap', cdict, 1024)
@np.vectorize
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
def Np(Array,number):
    unique, counts = np.unique(Array, return_counts=True)
    return dict(zip(unique, counts))[number]
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
def Distance(S1,S2):
    XS1,YS1 = np.mean(S1[0::2]),np.mean(S1[1::2])
    XS2,YS2 = np.mean(S2[0::2]),np.mean(S2[1::2])
    return np.max([np.sqrt((XS2-XS1)**2+(YS2-YS1)**2),10**(-5)])
def MakeAdjacency(data):
    Arr = np.zeros((data.shape[0],data.shape[0]),dtype=bool)
    #A = csr_matrix((data.shape[0],data.shape[0]))
    #A = np.array([[Distance(data[i],data[j]) for i in range(data.shape[0])] for j in range(data.shape[0])])
    #A = sorted( set( A.flatten() ) )[-2]/A
    for i in range(data.shape[0]):
        for j in range(data.shape[0]):
            if IsNeighbors(data[i],data[j]) or i==j:
                Arr[i,j] = True
    A = csr_matrix(Arr)
    return A
def MakeCluster(AdjacencyMatrix,resolution=0.2):
    louvain = skn.clustering.Louvain(resolution=resolution)
    labels = louvain.fit_transform(AdjacencyMatrix)
    #propagation = skn.clustering.PropagationClustering()
    #labels = propagation.fit_transform(AdjacencyMatrix)
    return labels
def VerifyClustering(Adjacency,labels,confidence=10):
    ddict = defaultdict(list)
    ToMerge = list()
    for i,j in zip(labels,np.arange(labels.shape[0])):
        ddict[i].append(j)
    for n,c1 in enumerate(ddict.keys()):
        for c2 in np.asarray(list(ddict.keys()))[n:]:
            if c1!=c2:
                Path = list()
                Nsample = min([confidence,len(ddict[c1]),len(ddict[c2])])
                #for _ in range(min([confidence,len(ddict[c1]),len(ddict[c2])])):
                #    i1,i2 = random.choice(ddict[c1]),random.choice(ddict[c2])
                #for i1 in ddict[c1]:
                #    for i2 in ddict[c2]:
                for i1 in random.sample(ddict[c1],Nsample):
                    for i2 in random.sample(ddict[c2],Nsample):
                        Path.append(set(skn.path.shortest_path(Adjacency,i1,i2)))
                if len(set.intersection(*Path)) == 0 and all([len(p)!=0 for p in Path]):
                    ToMerge.append((c1,c2))
    return ToMerge

def Merge(labels,ToMerge):
    for ij in ToMerge:
        labels[np.where(labels==ij[1])] = ij[0]
    #return labels
def GetAvSize(filename,Check=False,resolution = 0.2,Safe = False,confidence = 10):
    try :
        data = np.loadtxt(filename,dtype=float)
    except :
        print('no file' + filename)
        return (0,0)
    A = MakeAdjacency(data)
    labels = MakeCluster(A,resolution = resolution)
    if Safe :
        ToMerge = VerifyClustering(A,labels,confidence)
        Merge(labels,ToMerge)
    AvSize = np.mean([Np(labels,lab) for lab in set(labels)])
    if Check :
        for n,S in enumerate(data):
            X,Y = np.mean(S[0::2]),np.mean(S[1::2])
            plt.scatter(X,Y,color=cm((labels[n]-min(labels))/len(set(labels))))
        plt.show()
    return AvSize,len(set(labels))
