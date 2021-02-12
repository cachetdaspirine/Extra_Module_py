########################################################
## Functions which are specific to the model ###########
########################################################
import time
import numpy as np
from scipy.optimize import fmin_cg
import math
########################################################
########################################################

########################################################
########################################################
from ModelFunctions_v1 import *
########################################################
########################################################

########################################################
########################################################
def CalcualteTotalEnergy(xy_dof, *args):
    ## Calculates total energy from the values
    ## of degrees of freedom
    ################################################
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon, w_temp)= args
    ################################################    
    ## Frist determine the edge coordinates from DOF
    (xy_nodes_red_array_t, \
     xy_nodes_blue_array_t, tt_temp) = \
     XYfromDOF(xy_dof, w_temp)
    ################################################
    ##################
    N_triangles = 2*w_temp
    energy_temp = 0.0
    for ind_tr in range(N_triangles):
        level_down = math.floor(ind_tr/2)
        xyr = np.zeros([3, 2])
        xyb = np.zeros([3, 2])
        if ind_tr%2==0:
            ind_t1 = 2*level_down
            ind_t2 = (level_down%2) + 2*(level_down + 1)
            ind_t3 = 2*level_down + 1
        else:
            ind_t1 = 2*(level_down+1)
            ind_t2 = ((level_down+1)%2) + 2*(level_down) 
            ind_t3 = 2*(level_down+1) + 1
        #print(level_down, ind_tr, ind_t1, ind_t2, ind_t3)
        xyr[0, 0] = xy_nodes_red_array_t[ind_t1, 0]
        xyr[0, 1] = xy_nodes_red_array_t[ind_t1, 1]
        xyr[1, 0] = xy_nodes_red_array_t[ind_t2, 0]
        xyr[1, 1] = xy_nodes_red_array_t[ind_t2, 1]
        xyr[2, 0] = xy_nodes_red_array_t[ind_t3, 0]
        xyr[2, 1] = xy_nodes_red_array_t[ind_t3, 1]         
        xyb[0, 0] = xy_nodes_blue_array_t[ind_t1, 0]       
        xyb[0, 1] = xy_nodes_blue_array_t[ind_t1, 1]
        xyb[1, 0] = xy_nodes_blue_array_t[ind_t2, 0]
        xyb[1, 1] = xy_nodes_blue_array_t[ind_t2, 1]
        xyb[2, 0] = xy_nodes_blue_array_t[ind_t3, 0]
        xyb[2, 1] = xy_nodes_blue_array_t[ind_t3, 1]
        energy_temp = energy_temp + \
                      CalculateParticleEnergy(xyr, xyb, *args)
    energyTotal = energy_temp/float(N_triangles)
    #################################################
    return energyTotal
########################################################
########################################################

########################################################
########################################################
def XYfromDOF(xy_dof_t, w_t):
    ## Determine DOF values from edge coordinates
    ##
    Nnodes = int(2.0*float(w_t+1))
    ## Number of degree of freedom
    Nfree_level = int(math.floor((w_t+1)/2))
    Nfree = int(2*math.floor((w_t+1)/2) + 1)
    xy_nodes_red_array_temp = np.zeros([Nnodes, 2])
    xy_nodes_blue_array_temp = np.zeros([Nnodes, 2])
    tt_temp = xy_dof_t[Nfree-1]
    for ind_d in range(Nfree_level):
        ind_r = 2*ind_d
        ind_b = 2*ind_d + 1
        ##
        ind_rev = 2*(Nfree_level - 1 - ind_d)
        ind_1 = (ind_rev)
        ind_2 = ind_1 + 1
        ind_3 = Nnodes - 2 - ind_rev
        ind_4 = ind_3 + 1
        ##
        xy_nodes_red_array_temp[ind_1, 1] = - xy_dof_t[ind_r]
        xy_nodes_red_array_temp[ind_2, 1] = - xy_dof_t[ind_r]
        xy_nodes_red_array_temp[ind_3, 1] = xy_dof_t[ind_r]
        xy_nodes_red_array_temp[ind_4, 1] = xy_dof_t[ind_r]
        xy_nodes_red_array_temp[ind_1, 0] = float((ind_1/2)%2)*tt_temp/2.0
        xy_nodes_red_array_temp[ind_2, 0] = xy_nodes_red_array_temp[ind_1, 0] + \
                                            tt_temp
        xy_nodes_red_array_temp[ind_3, 0] = float((ind_3/2)%2)*tt_temp/2.0
        xy_nodes_red_array_temp[ind_4, 0] = xy_nodes_red_array_temp[ind_3, 0] + \
                                            tt_temp
        ##
        xy_nodes_blue_array_temp[ind_1, 1] = - xy_dof_t[ind_b]
        xy_nodes_blue_array_temp[ind_2, 1] = - xy_dof_t[ind_b]
        xy_nodes_blue_array_temp[ind_3, 1] = xy_dof_t[ind_b]
        xy_nodes_blue_array_temp[ind_4, 1] = xy_dof_t[ind_b]
        xy_nodes_blue_array_temp[ind_1, 0] = xy_nodes_red_array_temp[ind_1, 0]
        xy_nodes_blue_array_temp[ind_2, 0] = xy_nodes_red_array_temp[ind_2, 0]
        xy_nodes_blue_array_temp[ind_3, 0] = xy_nodes_red_array_temp[ind_3, 0]
        xy_nodes_blue_array_temp[ind_4, 0] = xy_nodes_red_array_temp[ind_4, 0]
        ##
    if w_t%2==0:
        ind_middle = int(w_t)
        xy_nodes_red_array_temp[ind_middle, 0] = float((ind_middle/2)%2)*tt_temp/2.0
        xy_nodes_red_array_temp[ind_middle + 1, 0] = xy_nodes_red_array_temp[ind_middle, 0] + \
                                                     tt_temp
        xy_nodes_blue_array_temp[ind_middle, 0] = xy_nodes_red_array_temp[ind_middle, 0]
        xy_nodes_blue_array_temp[ind_middle + 1, 0] = xy_nodes_red_array_temp[ind_middle + 1, 0]
    ##########
    return (xy_nodes_red_array_temp, \
            xy_nodes_blue_array_temp, tt_temp)
########################################################
########################################################

########################################################
########################################################
def InitialDOFlist(w_t, l0_red, l0_blue):
    ## Number of nodes
    Nnodes = int(2.0*float(w_t+1))
    ## Number of degree of freedom
    Nfree_level = int(math.floor((w_t+1)/2))
    Nfree = int(2*math.floor((w_t+1)/2) + 1)
    ## values of DOF
    xy_dof_at = np.zeros(Nfree)
    if w_t%2==0:
        y0_t_red = math.sqrt(3.0/4.0)*l0_red
        y0_t_blue = math.sqrt(3.0/4.0)*l0_blue
    elif w_t%2==1:        
        y0_t_red = math.sqrt(3.0/4.0)*l0_red/2.0
        y0_t_blue = math.sqrt(3.0/4.0)*l0_blue/2.0
    for ind_f in range(Nfree_level):
        xy_dof_at[2*ind_f] = y0_t_red+ \
                               ind_f*math.sqrt(3.0/4.0)*l0_red
        xy_dof_at[2*ind_f+1] = y0_t_blue+ \
                               ind_f*math.sqrt(3.0/4.0)*l0_blue
    xy_dof_at[Nfree-1] = (l0_red + l0_blue)/2.0
    ###########################
    return(xy_dof_at)
########################################################
########################################################

########################################################
########################################################
def FunctionFilamentEnergy(width, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    ###############
    args = (k_red, k_blue, k_soft, \
            k_Area_red, k_Area_blue, \
            l0_0, l0_red, l0_blue, l0_soft, \
            area_0_red, area_0_blue, \
            epsilon, width)
    ###############
    xy_dof_0 = InitialDOFlist(width, l0_red, l0_blue)
    ###############
    res = fmin_cg(CalcualteTotalEnergy, xy_dof_0, \
                  args=args, full_output=1, disp=0)
    res_min = res[0]
    ###############
    if res[4]!=0:
        print('\n')
        print('Error: There is an error in the filament energy calculations in the first sturcture')
        print('\n')
        print(res[4])
    ###############
    en_w_temp = CalcualteTotalEnergy(res_min, *args)
    ###############
    return en_w_temp
########################################################
########################################################

########################################################
########################################################
