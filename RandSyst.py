import os
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import time
import numpy as np
import pathlib
from ctypes import cdll
from ctypes import c_double
from ctypes  import c_long
from ctypes import c_int
from ctypes import c_bool
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_char_p
import copy

lib1 = cdll.LoadLibrary(
    str(pathlib.Path(__file__).parent.absolute()) + '/libRand.so')

lib1.CreateSystem.restype = POINTER(c_void_p)
lib1.CreateSystem.argtypes = [POINTER(c_int), POINTER(c_double), POINTER(c_double), c_int, c_int]
lib1.DeleteSystem.argtypes = [POINTER(c_void_p)]
lib1.CopySystem.argtypes = [POINTER(c_void_p)]
lib1.CopySystem.restype = POINTER(c_void_p)

lib1.UpdateSystemEnergy.argtypes = [
    POINTER(c_void_p), POINTER(c_int), c_int, c_int]
lib1.GetSystemEnergy.restype = c_double
lib1.GetSystemEnergy.argtypes = [POINTER(c_void_p)]

lib1.OutputSystemSite.argtypes = [POINTER(c_void_p), c_char_p]
lib1.OutputSystemSiteExtended.argtypes = [POINTER(c_void_p), c_char_p]

lib1.GetBulkEnergy.argtypes = [POINTER(c_void_p)]
lib1.GetBulkEnergy.restype = c_double


#lib1.AffineDeformation.argtypes = [POINTER(c_void_p),c_double,c_double]
#lib1.AffineDeformation.restype = c_double
lib1.AffineDeformation.argtypes = [POINTER(c_void_p),c_double,c_double,POINTER(c_int),c_int]
lib1.AffineDeformation.restype = c_double

lib1.Extension.argtypes = [POINTER(c_void_p),c_int]
lib1.Extension.restype = c_double

lib1.GetNDOF.argtypes = [POINTER(c_void_p)]
lib1.GetNDOF.restypes = c_int

lib1.GetHessian.argtypes = [POINTER(c_void_p), POINTER(c_double), c_int]
lib1.GetGradient.argtypes = [POINTER(c_void_p), POINTER(c_double), c_int]

lib1.GetDOFIndex.argtypes = [POINTER(c_void_p),POINTER(c_double)]#,POINTER(c_long),POINTER(c_long),POINTER(c_long),POINTER(c_long)]

lib2 = cdll.LoadLibrary(
    str(pathlib.Path(__file__).parent.absolute()) + '/libRand_expansion.so')

lib2.CreateSystem.restype = POINTER(c_void_p)
lib2.CreateSystem.argtypes = [POINTER(c_int), POINTER(c_double), POINTER(c_double), c_int, c_int]
lib2.DeleteSystem.argtypes = [POINTER(c_void_p)]
lib2.CopySystem.argtypes = [POINTER(c_void_p)]
lib2.CopySystem.restype = POINTER(c_void_p)

lib2.UpdateSystemEnergy.argtypes = [
    POINTER(c_void_p), POINTER(c_int), c_int, c_int]
lib2.GetSystemEnergy.restype = c_double
lib2.GetSystemEnergy.argtypes = [POINTER(c_void_p)]

lib2.OutputSystemSite.argtypes = [POINTER(c_void_p), c_char_p]

lib2.GetBulkEnergy.argtypes = [POINTER(c_void_p)]
lib2.GetBulkEnergy.restype = c_double


lib2.AffineDeformation.argtypes = [POINTER(c_void_p),c_double,c_double,POINTER(c_int),c_int]
lib2.AffineDeformation.restype = c_double

lib2.Extension.argtypes = [POINTER(c_void_p),c_int]
lib2.Extension.restype = c_double

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

class System:
    def __init__(self,CouplingMatrix=None, q0=None, State=None,old_system=None):
        if old_system == None:
                self.None_Copy(State, CouplingMatrix, q0)
        else:
            self.Copy(old_system)

    def None_Copy(self, State, CouplingMatrix, q0):
        # The system is created by a 2D array for a 2D system the shape[0] is
        # the  Y  lengft  while  the  shape[1]  is  the  X  shape  the   most
        # important part of this object is self.adress which is the adress of
        # the pointer toward the cpp object. Each time we call a c++ function
        # we have to give it the adress of the  pointer,  that  the  function
        # will interpret as a pointer toward the c++ object
        if q0.shape[0]==12:
            self.lib = lib2
        elif q0.shape[0]==9:
            self.lib = lib1
        self.Lx = State.shape[0]  # X size of the system !!!!!
        self.Ly = State.shape[1]  # Y size of the system !!!!!
        # --------------Convert the array into a pointer array---------------
        self.state = State  # store the value of the binary system as a 2D array
        array = np.zeros(State.shape[0] * State.shape[1], dtype=int)
        for i in range(State.shape[0]):
            for j in range(State.shape[1]):
                array[i + j * State.shape[0]] = State[i, j]
        # declare a pointer array of integer (ctypes type)
        Arraycpp = array.ctypes.data_as(POINTER(c_int))
        for i in range(array.shape[0]):
            # store all the array into this pointer array
            Arraycpp[i] = array[i]
        # -------------------------------------------------------------------
        # store the value of the elastic parameters
        # -------------------------------------------------------------------
        self.CouplingMatrix = CouplingMatrix
        self.q0 = q0

        MC = copy.copy(self.CouplingMatrix.flatten())
        MCcpp = MC.ctypes.data_as(POINTER(c_double))
        for i in range(MC.shape[0]):
            MCcpp[i] = MC[i]

        Q0 = copy.copy(self.q0)
        Q0cpp = Q0.ctypes.data_as(POINTER(c_double))
        for i in range(Q0.shape[0]):
            Q0cpp[i] = Q0[i]

        self.ActualizeNp()  # keep track of the number of particle (number of 1) in the system
        # ---------------------Create the cpp object-------------------------
        # create the system, all the argument are require here !!!!
        self.Adress = self.lib.CreateSystem(
            Arraycpp, MCcpp,Q0cpp,self.Lx,self.Ly)
        # --------------------Store the value of the Energy------------------
        # store the value of the Energy (get energy only returns a number and doesn't reactualize the equilibrium of the system).
        self.Energy = self.lib.GetSystemEnergy(self.Adress)
        self.MAP = cm


    def Copy(self, old_system):
        self.Lx = old_system.Lx
        self.Ly = old_system.Ly
        self.state = old_system.state
        self.CouplingMatrix = old_system.CouplingMatrix
        self.q0 = old_system.q0
        if self.q0==12:
            self.lib = lib2
        elif self.q0==9:
            self.lib = lib1
        self.ActualizeNp()
        self.Adress = self.lib.CopySystem(old_system.Adress)
        self.Energy = self.lib.GetSystemEnergy(self.Adress)

    def __del__(self):
        # deleting pointers is important in c++
        #print('delete the system')
        self.lib.DeleteSystem(self.Adress)

    def get_Hessian(self):
        NDOF = self.lib.GetNDOF(self.Adress)
        Hessian = np.zeros(NDOF*NDOF,dtype=np.double)
        #Arraycpp =
        self.lib.GetHessian(self.Adress, Hessian.ctypes.data_as(POINTER(c_double)) , NDOF)
        return np.reshape(Hessian,(-1,NDOF))
    def PlotHessVect(self,Vect):
        Hess = self.get_Hessian()
        Mapping = self.get_IndexList()
        self.PrintPerSite('ToPlot.txt')
        Data = np.loadtxt('ToPlot.txt', dtype=float)
        os.system('rm -rf ToPlot.txt')

        plt.quiver()
    def get_IndexList(self):
        NDOF = self.lib.GetNDOF(self.Adress)
        #Is = np.zeros(NDOF,dtype=np.int)
        #Js = np.zeros(NDOF,dtype=np.int)
        #Ks = np.zeros(NDOF,dtype=np.int)
        #Xs = np.zeros(NDOF,dtype=np.int)
        XYs = np.zeros(NDOF,dtype=np.double)
        self.lib.GetDOFIndex(self.Adress,XYs.ctypes.data_as(POINTER(c_double)))
                            #Is.ctypes.data_as(POINTER(c_long)),
                            #Js.ctypes.data_as(POINTER(c_long)),
                            #Ks.ctypes.data_as(POINTER(c_long)),
                            #Xs.ctypes.data_as(POINTER(c_long)))
        #Maps = np.loadtxt('Mapping.txt',dtype=int)
        #os.system('rm Mapping.txt')
        return XYs

    def get_Gradient(self):
        NDOF = self.lib.GetNDOF(self.Adress)
        length= 2 * NDOF
        Gradient = np.zeros(length,dtype=np.double)
        self.lib.GetGradient(self.Adress,Gradient.ctypes.data_as(POINTER(c_double)),length)
        return Gradient

    def GetBulkEnergy(self):
        return self.lib.GetBulkEnergy(self.Adress)
    def Evolv(self, NewState):
        self.ActualizeNp()
        # ------------Convert the new state into a pointer array-------------
        if NewState:
            self.state = NewState
        else :
            NewState=self.state
        array = np.zeros(self.state.shape[0] * self.state.shape[1], dtype=int)
        for i in range(self.state.shape[0]):
            for j in range(self.state.shape[1]):
                array[i + j * self.state.shape[0]] = self.state[i, j]
        Arraycpp = array.ctypes.data_as(POINTER(c_int))
        for i in range(array.shape[0]):
            Arraycpp[i] = array[i]
        # -------------------------------------------------------------------
        # Second most important function you give a new state and it computes
        # the new equilibrium  state,  ît just goes faster as the newstate is
        # close to the previous one.
        if NewState.shape[0] != self.Lx or NewState.shape[1] != self.Ly:
            # if we changed the size of the system, we remake the whole system
            self.Lx = NewState.shape[0]
            self.Ly = NewState.shape[ju1]
            self.lib.DeleteSystem(self.Adress)

            MC = copy.copy(self.CouplingMatrix.flatten())
            MCcpp = MC.ctypes.data_as(POINTER(c_double))
            for i in range(MC.shape[0]):
                MCcpp[i] = MC[i]

            Q0 = copy.copy(self.q0)
            Q0cpp = Q0.ctypes.data_as(POINTER(c_double))
            for i in range(Q0.shape[0]):
                Q0cpp[i] = Q0[i]

            self.Adress = self.lib.CreateSystem(
                Arraycpp, MCcpp,Q0cpp,self.Lx, self.Ly)
            self.Energy = self.lib.GetSystemEnergy(self.Adress)
            print('create a new system')
        else:
            self.lib.UpdateSystemEnergy(
                self.Adress, Arraycpp, self.Lx, self.Ly)
            self.Energy = self.lib.GetSystemEnergy(self.Adress)

    def PrintPerSite(self, Name='NoName.txt',Extended = False):
        # output the sytem per site (easier if you wanna plot the sites).
        if self.Np < 1:
            print("can t output an empty system")
            return 0.
        if Extended :
            self.lib.OutputSystemSiteExtended(self.Adress, Name.encode('utf-8'))
        else:
            self.lib.OutputSystemSite(self.Adress, Name.encode('utf-8'))
    def get_sites_position_as_array(self):
        self.PrintPerSite('to_plot.txt',Extended=True)
        Data = np.loadtxt('to_plot.txt')
        os.system('rm to_plot.txt')
        array_of_position = np.zeros(self.state.shape,dtype=np.ndarray)
        for site in Data:
            array_of_position[int(site[-2]),int(site[-1])] = site[:-3]
        return array_of_position
    def GetNodePerSite(self):
        self.PrintPerSite('ToReturn.txt')
        Data = np.loadtxt('ToReturn.txt',dtype=float)
        os.system('rm -f ToReturn.txt')
        return Data

    def PlotPerSite(self, figuresize=(7, 5), Zoom=1.,Fill=True,FillColor='stress',FIGAX=None,ec=(0,0,0,1)):
        # this one has a trick, it only 'works' on UNIX system and
        # it requires to be autorized to edit and delete file. The
        # idea is to use the function  in  order  to  PrintPersite
        # create  a  file  that we load, then delete. Then  we use
        # matplotlib triangle patches to plot the system
        if self.Np < 1:
            print("can t output an empty system")
            return 0.
        # Directly plot the whole system as patches of polygon, it just require to save a file that it will delete
        if FIGAX:
            fig,ax = FIGAX
        else:
            fig, ax = plt.subplots(figsize=figuresize)
        self.PrintPerSite('ToPlot.txt')
        Data = np.loadtxt('ToPlot.txt', dtype=float)
        os.system('rm -rf ToPlot.txt')
        XC, YC = 0, 0
        if type(Data[0]) != np.ndarray:
            Data = np.array([Data])
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
        C = (np.array(C)/max(C))#-0.5)*2
        if FillColor == 'plain':
            for n,XY in enumerate(Hex):
                ax.add_patch(Polygon(XY, closed=True, linewidth=0.8, fill=Fill, fc=(
                    0.41, 0.83, 0.94, 0.5), ec=(0, 0, 0, 1), ls='-', zorder=0))
        else:# C = 'stress':
            if Fill:
                coll = PolyCollection(Hex,array=C, closed=True, linewidth=0.8, cmap=cm, ec=(0,0,0,1), ls='-', zorder=0)
            else :
                coll = PolyCollection(Hex,facecolors = 'none', closed=True, linewidth=0.8, cmap=cm, ec=(0,0,0,1), ls='-', zorder=0)
            ax.add_collection(coll)
            fig.colorbar(coll, ax=ax)
            #ax.add_patch(Polygon(XY, closed=True, linewidth=0.8, fill=Fill, fc=cm(C[n]), ec=ec, ls='-', zorder=0))
        ax.set_aspect(aspect=1.)
        ax.set_xlim([XC / Data.shape[0] - 1 / Zoom * np.sqrt(Data.shape[0]),
                     XC / Data.shape[0] + 1 / Zoom * np.sqrt(Data.shape[0])])
        ax.set_ylim([YC / Data.shape[0] - 1 / Zoom * np.sqrt(Data.shape[0]),
                     YC / Data.shape[0] + 1 / Zoom * np.sqrt(Data.shape[0])])

        # plt.show()
        return fig, ax

    def ActualizeNp(self):
        # transform the array of 0 and  1 into a dictionnary, which key is
        # the 0s or 1s and the respective value is the number of particles
        # or the number of empty sites
        try:
            unique, counts = np.unique(self.state, return_counts=True)
            self.Np = dict(zip(unique, counts))[1]
        except:
            self.Np = 0

    def Extension(self,ax):
        if any([s==0 for S in self.state for s in S]):
            print('the system isn t full of particle, we can t apply any deformation')
            raise ValueError
        return self.lib.Extension(self.Adress,ax)

    def AffineDeformation(self, deformation_x = 0, deformation_y = 0, NodesIndex = [0,1]):
        if type(NodesIndex) != np.ndarray:
            NodesIndex= np.array(NodesIndex)
        Object=NodesIndex.ctypes.data_as(POINTER(c_int))
        for i in range(NodesIndex.shape[0]):
            Object[i] = NodesIndex[i]
        # this is not compatible with the spring expansion model
        return self.lib.AffineDeformation(self.Adress,deformation_x,deformation_y,
                                          Object,
                                          NodesIndex.shape[0])
