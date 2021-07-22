#!/home/hugo/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
import matplotlib.lines as mlines
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
import sys
from matplotlib.colors import LinearSegmentedColormap

import math
def PlotTriangle(*argv):
    print(argv)
    def truncate(number, digits) -> float:
        stepper = 10.0 ** digits
        return math.trunc(stepper * number) / stepper

    if len(argv)<1 :
        print("File name not specified")
        sys.exit()

    if type(argv[0]) != str :
        print("File name not specified")
        sys.exit()

    # Start by sorting the argument into several curves:
    for n,i in enumerate(argv):
        if n==0:
            argument=list()
            continue
        else:
            argument.append(i)

    Edge=True
    zoom=1
    for n,c in enumerate(argument):
        if c=='NoEdge':
            Edge=False
        if c=='zoom':
            zoom=float(argument[n+1])


    cdict = {'blue':   ((0.0,  0.9,0.9),
                        (0.5,  0.4, 0.4),
                        (1.0,  0.1, 0.1)),

             'green': ((0.0,  0.5, 0.5),
                       (0.5 , 1, 1),
                       (1.0,  0.3, 0.3)),

             'alpha': ((0.0,  1, 1),
                       (0.5 , 0.8, 0.8),
                       (1.0,  1, 1)),

             'red':  ((0.0,  0.4, 0.4),
                       (0.5,  0.5, 0.5),
                       (1.0,  0.9,0.9)),
    }
    cm = LinearSegmentedColormap('my_colormap', cdict, 1024)
    fig,ax=plt.subplots()#figsize=(7,5))

    Data=np.loadtxt(argv[0],dtype=np.float64)

    if Data.shape[1] == 6:
        Type='Triangle'
        Stress=False
    elif len(Data[1]) == 12:
        Type='Hexagon'
        Stress = False
    elif len(Data[1]) == 13:
        Type='Hexagon'
        Stress= True
    else :
        print('unknown particle type')

    XC,YC=0,0

    # Start to plot the triangle themself
    if Stress:
        Hex=list()
        C = list()
        for ligne in Data:
            XY = []
            for i in range((ligne.shape[0]-1) // 2):
                XY.append([ligne[2 * i], ligne[2 * i + 1]])
            XC += sum(np.transpose(XY)[0]) / len(XY)
            YC += sum(np.transpose(XY)[1]) / len(XY)
            Color = ligne[-1]
            Hex.append(XY)
            C.append(Color)
        C = (np.array(C)/max(C)-0.5)*2
        #for n,XY in enumerate(Hex):
        #   ax.add_patch(Polygon(XY, closed=True, linewidth=0.8, fill=True, fc=cm(C[n]), ec=(0,0,0,1), ls='-', zorder=0))
        coll = PolyCollection(Hex,array=C, closed=True, linewidth=0.8, cmap=cm, ec=(0,0,0,1), ls='-', zorder=0)
        ax.add_collection(coll)
        #cax = fig.add_axes([0,0,1,1])
        fig.colorbar(coll, ax=ax)
    else :
        for ligne in Data :
            XY=[]
            for i in range((ligne.shape[0]+1)//2):
                XY.append([ligne[2*i],ligne[2*i+1]])
            XC+=sum(np.transpose(XY)[0])/len(XY)
            YC+=sum(np.transpose(XY)[1])/len(XY)
            ax.add_patch(Polygon(XY,closed=True,linewidth=0.8,fill=True,fc=(0.41,0.83,0.94,0.5),ec=(0,0,0,1),ls='-',zorder=0))#fc=(0.,0.,0.,0.),ec=(0,0,0,1),ls='-',zorder=0))
    # Plot the free edge in red
    if Edge:
        edges=list()
        #----------Make a list with all the edge in my system---------
        for ligne in Data:
            PointsX=[truncate(ligne[2*i],3) for i in range(ligne.__len__()//2)]
            PointsY=[truncate(ligne[2*i+1],3) for i in range(ligne.__len__()//2)]
            Points=list(zip(PointsX,PointsY))
            for i in range(Points.__len__()):
                if Type=='Triangle':
                    edges.append(sorted([Points[i],Points[(i+1)%3]]))
                elif Type=='Hexagon':
                    edges.append(sorted([Points[i],Points[(i+1)%6]]))
        #----------Convert this list into a dictionnary which value is the number of appearance
        unique, counts = np.unique(edges, return_counts=True,axis=0)
        FreeEdgeIndex=np.where(counts==1)
        for n,index in enumerate(FreeEdgeIndex[0]):
            X , Y =list(),list()
            for Points in unique[index]:
                X.append(Points[0])
                Y.append(Points[1])
            ax.add_line(mlines.Line2D(X,Y,color=(231./255.,81./255.,19./255.)))

    ax.set_aspect(aspect=1.)
    ax.set_xlim([XC/Data.shape[0]-np.sqrt(Data.shape[0])/zoom,XC/Data.shape[0]+np.sqrt(Data.shape[0])/zoom])
    ax.set_ylim([YC/Data.shape[0]-np.sqrt(Data.shape[0])/zoom,YC/Data.shape[0]+np.sqrt(Data.shape[0])/zoom])

    return fig,ax
