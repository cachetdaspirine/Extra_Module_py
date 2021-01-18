#!/home/hugo/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines
from matplotlib.colors import LinearSegmentedColormap
import Instruction as In

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
import math
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
cm = LinearSegmentedColormap('my_colormap', cdict, 1024)

def Plot(File,Figsize=(7,5),zoom=1,Edge=True,axe=False):
    if axe:
        ax=axe
        fig=plt.figure(figsize=Figsize)
    else:
        fig,ax=plt.subplots(figsize=Figsize)

    Data=np.loadtxt(File,dtype=float)

    XC,YC=0,0

    for ligne in Data :
        XY=[]
        for i in range(ligne.shape[0]//2):
            XY.append([ligne[2*i],ligne[2*i+1]])
        XC+=sum(np.transpose(XY)[0])/len(XY)
        YC+=sum(np.transpose(XY)[1])/len(XY)
        ax.add_patch(Polygon(XY,closed=True,linewidth=0.8,fill=True,fc=(0.41,0.83,0.94,0.5),ec=(0,0,0,1),ls='-',zorder=0))
    if Edge:
        edges=list()
        #----------Make a list with all the edge in my system---------
        for ligne in Data:
            PointsX=[truncate(ligne[2*i],3) for i in range(ligne.__len__()//2)]
            PointsY=[truncate(ligne[2*i+1],3) for i in range(ligne.__len__()//2)]
            Points=list(zip(PointsX,PointsY))
            for i in range(Points.__len__()):
                edges.append(sorted([Points[i],Points[(i+1)%3]]))
        #----------Convert this list into a dictionnary which value is the number of appearance
        unique, counts = np.unique(edges, return_counts=True,axis=0)
        FreeEdgeIndex=np.where(counts==1)
        for n,index in enumerate(FreeEdgeIndex[0]):
            X , Y =list(),list()
            for Points in unique[index]:
                X.append(Points[0])
                Y.append(Points[1])
            ax.add_line(mlines.Line2D(X,Y,color=(231./255.,81./255.,19./255.)))

    ax.set_xlim([XC/Data.shape[0]-(2*np.sqrt(Data.shape[0]))/zoom,XC/Data.shape[0]+(2*np.sqrt(Data.shape[0]))/zoom])
    ax.set_ylim([YC/Data.shape[0]-(2*np.sqrt(Data.shape[0]))/zoom,YC/Data.shape[0]+(2*np.sqrt(Data.shape[0]))/zoom])

    return fig,ax
