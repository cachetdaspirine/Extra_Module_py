import numpy as np
import Aggregate as agg

def MakeCluster(Array):
    #make a list of list of object aggregate.
    #create a list of particle to make an aggregate of it.
    Pos=Array
    NpartTot=1
    Aggregates=np.array([])
    while Pos.shape[0]!=0:
        #print(Pos.shape[0])
        #Position1Agg=np.array([Pos[0,:]])
        Position1Agg=np.array([])
        Tampon=np.array([Pos[0,:]])
        Pos=np.delete(Pos,0,axis=0)
        while Tampon.shape[0]!=0:
            #put all the node of the tampon in an array
            node=np.array([])
            for i in range(Tampon.shape[0]):
                if node.shape[0]!=0 :
                    node=np.append(node,[Tampon[i,3:5],Tampon[i,9:11],Tampon[i,15:17]],axis=0)
                else :
                    node =np.array([Tampon[i,3:5],Tampon[i,9:11],Tampon[i,15:17]])

            TamponOldLength=Tampon.shape[0]
            for i in range(node.shape[0]):
                j=0
                while j<Pos.shape[0]:
                    if all(Pos[j,3:5]==node[i]) or all(Pos[j,9:11]==node[i]) or all(Pos[j,15:17]==node[i]):
                        Tampon=np.append(Tampon,[Pos[j,:]],axis=0)
                        Pos=np.delete(Pos,j,axis=0)
                        NpartTot+=1
                    else:
                        j+=1

            for i in range(TamponOldLength):
                if Position1Agg.shape[0]!=0 :
                    Position1Agg=np.append(Position1Agg,[Tampon[0,:]],axis=0)
                else :
                    Position1Agg=np.array([Tampon[0,:]])
                Tampon=np.delete(Tampon,0,0)

        #print(Aggregates.shape[0])
        if Aggregates.shape[0]!=0:
            Aggregates=np.append(Aggregates,agg.Aggregate(Array=Position1Agg))
        else:
            Aggregates=np.array([agg.Aggregate(Array=Position1Agg)])
    return Aggregates
