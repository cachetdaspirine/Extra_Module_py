########################################################
## Filament, energy per particle is determined
## by conugate gradiant algorithm.
########################################################
import numpy as np
from scipy.optimize import fmin_cg
import math
########################################################
########################################################

########################################################
########################################################
########################################################
########################################################
## main program ########################################
########################################################

########################################################
### Other Functions ####################################
#from ModelFunctions_v1 import *
from FilamentsStructure_1_v2 import *
#from FilamentsStructure_2_v2 import *
########################################################
########################################################

########################################################
########################################################
def DetermineInfiniFilamentEnergy(M_ij, rho0_vec, J_surface, width_max):
    ########################################################
    ## Red is the larger triangle ##########################
    ########################################################
    ########################################################
    ########################################################
    args = (M_ij, rho0_vec)
    ########################################################
    ### First Determine the bulk energy ####################
    BulkE = CalculateBulkEnergy(*args)
    ########################################################
    ### First Determine the bulk energy ####################
    ## Determine the filament energies for given parameters
    width_list = np.zeros(width_max-1)
    energy_1_list = np.zeros(width_max-1)
    #energy_2_list = np.zeros(width_max-1)
    for ind_w in range(1, width_max):
        #print(ind_w)
        width = ind_w
        width_list[ind_w-1] = ind_w
        energy_1_list[ind_w-1] = FunctionFilamentEnergy_FirstConfiguation(width, *args)
        #energy_2_list[ind_w-1] = FunctionFilamentEnergy_SecondConfiguation(width, *args)
        ###############################
        ## Add surface energy per particle
        energy_1_list[ind_w-1] = energy_1_list[ind_w-1] + 4.0*J_surface/width
        #energy_2_list[ind_w-1] = energy_2_list[ind_w-1] + 4.0*J_surface/width
    ########################################################

        #print(ind_w, energy_1_list[ind_w-1])#, energy_2_list[ind_w-1])
        ###
    return (width_list, energy_1_list, BulkE)#, energy_2_list, BulkE)
########################################################
########################################################
