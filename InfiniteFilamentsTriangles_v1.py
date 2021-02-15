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
from FilamentsStructure_v1 import *
########################################################
########################################################

########################################################
########################################################
def DetermineInfiniFilamentEnergy(k_red, k_blue, k_soft, \
                                  k_Area_red, k_Area_blue, \
                                  l0_0, \
                                  epsilon, J_surface, width_max):
    ########################################################
    ## Red is the larger triangle ##########################
    ########################################################
    ########################################################
    l0_red = l0_0*(1.0 + epsilon)
    l0_blue = l0_0*(1.0 - epsilon)
    l0_soft = math.sqrt(1.0 + epsilon*epsilon/3.0)*l0_0
    area_0_red = math.sqrt(3.0)*l0_red*l0_red/4.0
    area_0_blue = math.sqrt(3.0)*l0_blue*l0_blue/4.0
    ########################################################
    args = (k_red, k_blue, k_soft, \
            k_Area_red, k_Area_blue, \
            l0_0, l0_red, l0_blue, l0_soft, \
            area_0_red, area_0_blue, \
            epsilon)
    ########################################################
    ### First Determine the bulk energy ####################
    BulkE = CalculateBulkEnergy(*args)
    ########################################################
    ### First Determine the bulk energy ####################
    ## Determine the filament energies for given parameters
    width_list = np.zeros(width_max-1)
    energy_list = np.zeros(width_max-1)
    for ind_w in range(1, width_max):
        #print(ind_w)
        width = ind_w
        width_list[ind_w-1] = ind_w
        energy_list[ind_w-1] = FunctionFilamentEnergy(width, *args)
        ###############################
        ## Add surface energy per particle
        energy_list[ind_w-1] = energy_list[ind_w-1] + J_surface/width
    ########################################################
    return (width_list, energy_list, BulkE)
########################################################
########################################################
