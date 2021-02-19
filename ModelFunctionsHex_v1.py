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
def CalculateParticleEnergy(edge_coordinate, *args):
    ###################################################
    ## Calcualtes the energy of a individual particles
    ## from coordinates of it edges.
    ## Note that edge numbered with "0" must be red.
    ###################################################
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    energy_sum_temp = 0.0
    ## spring energy
    energy_sum_spring_temp = 0.0
    for ind_sp in range(6):
        length_x = edge_coordinate[ind_sp, 0] - edge_coordinate[(ind_sp+1)%6, 0]
        length_y = edge_coordinate[ind_sp, 1] - edge_coordinate[(ind_sp+1)%6, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = 0.5*k_soft*(lengthSQ - l0_soft)*(lengthSQ - l0_soft)
        energy_sum_spring_temp = energy_sum_spring_temp + en_temptemp
        if (ind_sp)%2==0:
            k_color = k_red
            l0_color = l0_red
        elif (ind_sp)%2==1:
            k_color = k_blue
            l0_color = l0_blue
        length_x = edge_coordinate[ind_sp, 0] - edge_coordinate[(ind_sp+2)%6, 0]
        length_y = edge_coordinate[ind_sp, 1] - edge_coordinate[(ind_sp+2)%6, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = 0.5*k_color*(lengthSQ - l0_color)*(lengthSQ - l0_color)
        energy_sum_spring_temp = energy_sum_spring_temp + en_temptemp        
    ## Area term in energy
    ## Red Triangle
    area_red = 0.0
    for ind_ar in range(3):
        ind_ai = 2*ind_ar
        x1_temp = edge_coordinate[ind_ai, 0]
        y1_temp = edge_coordinate[ind_ai, 1]
        x2_temp = edge_coordinate[(ind_ai+2)%6, 0]
        y2_temp = edge_coordinate[(ind_ai+2)%6, 1]
        sum_temp_area = x1_temp*y2_temp - x2_temp*y1_temp
        area_red = area_red + 0.5*sum_temp_area
    if area_red<0:
        area_red = 1000000.0
    energy_sum_temp = energy_sum_spring_temp 
    ## Blue Triangle
    area_blue = 0.0
    for ind_ar in range(3):
        ind_ai = 2*ind_ar + 1
        x1_temp = edge_coordinate[ind_ai, 0]
        y1_temp = edge_coordinate[ind_ai, 1]
        x2_temp = edge_coordinate[(ind_ai+2)%6, 0]
        y2_temp = edge_coordinate[(ind_ai+2)%6, 1]
        sum_temp_area = x1_temp*y2_temp - x2_temp*y1_temp
        area_blue = area_blue + 0.5*sum_temp_area
    if area_blue<0:
        area_blue = 1000000.0
    energy_sum_temp = energy_sum_spring_temp \
                      + 0.5*(k_Area_red/area_0_red)*\
                      (area_red - area_0_red)*\
                      (area_red - area_0_red)\
                      + 0.5*(k_Area_blue/area_0_blue)*\
                      (area_blue - area_0_blue)*\
                      (area_blue - area_0_blue)
    #########
    return energy_sum_temp
########################################################
########################################################

########################################################
########################################################
def ForcesUnit(edge_coordinate, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############    
    force_unit_temp = np.zeros([6, 2])
    ## Area energy
    area_red_temp = 0.0
    area_blue_temp = 0.0
    for ind_ar in range(3):
        ##
        ind_ai = 2*ind_ar
        x1_temp = edge_coordinate[ind_ai, 0]
        y1_temp = edge_coordinate[ind_ai, 1]
        x2_temp = edge_coordinate[(ind_ai+2)%6, 0]
        y2_temp = edge_coordinate[(ind_ai+2)%6, 1]
        sum_temp_area = x1_temp*y2_temp - x2_temp*y1_temp
        area_red_temp = area_red_temp + 0.5*sum_temp_area
        ##
        ind_ai = 2*ind_ar + 1
        x1_temp = edge_coordinate[ind_ai, 0]
        y1_temp = edge_coordinate[ind_ai, 1]
        x2_temp = edge_coordinate[(ind_ai+2)%6, 0]
        y2_temp = edge_coordinate[(ind_ai+2)%6, 1]
        sum_temp_area = x1_temp*y2_temp - x2_temp*y1_temp
        area_blue_temp = area_blue_temp + 0.5*sum_temp_area
    ###########    
    for ind_fu in range(6):
        x1_temp = edge_coordinate[ind_fu, 0]
        y1_temp = edge_coordinate[ind_fu, 1]
        x2_temp = edge_coordinate[(ind_fu+1)%6, 0]
        y2_temp = edge_coordinate[(ind_fu+1)%6, 1]
        x0_temp = edge_coordinate[(ind_fu-1)%6, 0]
        y0_temp = edge_coordinate[(ind_fu-1)%6, 1]
        x3_temp = edge_coordinate[(ind_fu+2)%6, 0]
        y3_temp = edge_coordinate[(ind_fu+2)%6, 1]
        x4_temp = edge_coordinate[(ind_fu-2)%6, 0]
        y4_temp = edge_coordinate[(ind_fu-2)%6, 1]
        if (ind_fu)%2==0:
            k_color = k_red
            l0_color = l0_red
        elif (ind_fu)%2==1:
            k_color = k_blue
            l0_color = l0_blue
        length_s1 = math.sqrt((x1_temp - x2_temp)*(x1_temp - x2_temp) +\
                           (y1_temp - y2_temp)*(y1_temp - y2_temp))
        length_s2 = math.sqrt((x1_temp - x0_temp)*(x1_temp - x0_temp) +\
                           (y1_temp - y0_temp)*(y1_temp - y0_temp))
        length_c1 = math.sqrt((x1_temp - x3_temp)*(x1_temp - x3_temp) +\
                           (y1_temp - y3_temp)*(y1_temp - y3_temp))
        length_c2 = math.sqrt((x1_temp - x4_temp)*(x1_temp - x4_temp) +\
                           (y1_temp - y4_temp)*(y1_temp - y4_temp))
        force_unit_temp[ind_fu, 0] = k_soft*(x1_temp - x2_temp)*\
                                     (length_s1 - l0_soft)/length_s1+\
                                     k_soft*(x1_temp - x0_temp)*\
                                     (length_s2 - l0_soft)/length_s2+\
                                     k_color*(x1_temp - x3_temp)*\
                                     (length_c1 - l0_color)/length_c1 +\
                                     k_color*(x1_temp - x4_temp)*\
                                     (length_c2 - l0_color)/length_c2
        force_unit_temp[ind_fu, 1] = k_soft*(y1_temp - y2_temp)*\
                                     (length_s1 - l0_soft)/length_s1 +\
                                     k_soft*(y1_temp - y0_temp)*\
                                     (length_s2 - l0_soft)/length_s2+\
                                     k_color*(y1_temp - y3_temp)*\
                                     (length_c1 - l0_color)/length_c1 +\
                                     k_color*(y1_temp - y4_temp)*\
                                     (length_c2 - l0_color)/length_c2
        if ind_fu%2==0:
            force_unit_temp[ind_fu, 0] = force_unit_temp[ind_fu, 0] +\
                                         0.5*(k_Area_red/area_0_red)*\
                                         (y3_temp - y4_temp)*\
                                         (area_red_temp - area_0_red)
            force_unit_temp[ind_fu, 1] = force_unit_temp[ind_fu, 1] +\
                                         0.5*(k_Area_red/area_0_red)*\
                                         (x4_temp - x3_temp)*\
                                         (area_red_temp - area_0_red)
        elif ind_fu%2==1:
            force_unit_temp[ind_fu, 0] = force_unit_temp[ind_fu, 0] +\
                                         0.5*(k_Area_blue/area_0_blue)*\
                                         (y3_temp - y4_temp)*\
                                         (area_blue_temp - area_0_blue)
            force_unit_temp[ind_fu, 1] = force_unit_temp[ind_fu, 1] +\
                                         0.5*(k_Area_blue/area_0_blue)*\
                                         (x4_temp - x3_temp)*\
                                         (area_blue_temp - area_0_blue)
    ##        
    return (force_unit_temp)
########################################################
########################################################

########################################################
########################################################
def FunctionBulkEnergy(l_temp, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    l_ini = l_temp
    l_ini_60 = l_ini*math.sqrt(3.0)/2.0
    edge_6 = np.zeros([6,2])
    edge_6[0, 0] = 0.0
    edge_6[0, 1] = 0.0
    edge_6[1, 0] = l_ini
    edge_6[1, 1] = 0.0
    edge_6[2, 0] = l_ini*3.0/2.0
    edge_6[2, 1] = l_ini_60
    edge_6[3, 0] = l_ini
    edge_6[3, 1] = 2.0*l_ini_60
    edge_6[4, 0] = 0.0
    edge_6[4, 1] = 2.0*l_ini_60
    edge_6[5, 0] = -l_ini/2.0
    edge_6[5, 1] = l_ini_60
    ##
    eBulk_temp = CalculateParticleEnergy(edge_6, *args)
    ##
    return eBulk_temp
########################################################
########################################################

########################################################
########################################################
def GradFunctionBulkEnergy(l_temp, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    l_ini = l_temp
    l_ini_60 = l_ini*math.sqrt(3.0)/2.0
    ##
    edge_6 = np.zeros([6,2])
    edge_6[0, 0] = 0.0
    edge_6[0, 1] = 0.0
    edge_6[1, 0] = l_ini
    edge_6[1, 1] = 0.0
    edge_6[2, 0] = l_ini*3.0/2.0
    edge_6[2, 1] = l_ini_60
    edge_6[3, 0] = l_ini
    edge_6[3, 1] = 2.0*l_ini_60
    edge_6[4, 0] = 0.0
    edge_6[4, 1] = 2.0*l_ini_60
    edge_6[5, 0] = -l_ini/2.0
    edge_6[5, 1] = l_ini_60
    grad_bulk_list = ForcesUnit(edge_6, *args)
    ##
    grad_l_ini_60 = math.sqrt(3.0)/2.0
    grad_6 = np.zeros([6,2])
    grad_6[0, 0] = 0.0
    grad_6[0, 1] = 0.0
    grad_6[1, 0] = 1
    grad_6[1, 1] = 0.0
    grad_6[2, 0] = 3.0/2.0
    grad_6[2, 1] = grad_l_ini_60
    grad_6[3, 0] = 1.0
    grad_6[3, 1] = 2.0*grad_l_ini_60
    grad_6[4, 0] = 0.0
    grad_6[4, 1] = 2.0*grad_l_ini_60
    grad_6[5, 0] = -1.0/2.0
    grad_6[5, 1] = grad_l_ini_60
    ##
    grad_edge = 0.0
    for ind_g in range(6):
        grad_edge =  grad_edge +\
                    grad_bulk_list[ind_g, 0]\
                    *grad_6[ind_g, 0]
        grad_edge =  grad_edge +\
                    grad_bulk_list[ind_g, 1]\
                    *grad_6[ind_g, 1]
    #######
    return grad_edge
########################################################
########################################################

########################################################
########################################################
def CalculateBulkEnergy(*args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    res = fmin_cg(FunctionBulkEnergy, l0_soft, \
                  GradFunctionBulkEnergy, \
                  args=args, full_output=1, disp=0)
    res_m = res[0]
    if res[4]!=0:
        print('\n')
        print('Error: There is an error in the bulk energy calculations')
        print('\n')
    eBulk = FunctionBulkEnergy(res_m, *args)
    ###############
    return eBulk
########################################################
########################################################

########################################################
########################################################
