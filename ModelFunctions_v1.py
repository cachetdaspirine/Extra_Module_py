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
## This function finds an area of a triangle using
## the coordinates of its corners.
def FindAreaOfTriangles(coords_corners_1):
    area_triangle_temp_1 = 0.0
    ###############################
    x0_temp = coords_corners_1[0, 0]
    y0_temp = coords_corners_1[0, 1]
    x1_temp = coords_corners_1[1, 0]
    y1_temp = coords_corners_1[1, 1]
    x2_temp = coords_corners_1[2, 0]
    y2_temp = coords_corners_1[2, 1]
    area_triangle_temp_1 = 0.5*(x1_temp*y0_temp - x0_temp*y1_temp) + \
                           0.5*(x2_temp*y1_temp - x1_temp*y2_temp) + \
                           0.5*(x0_temp*y2_temp - x2_temp*y0_temp)
    return (abs(area_triangle_temp_1))
########################################################
########################################################

########################################################
########################################################
def CalculateParticleEnergy(edge_coordinate_1, \
                            edge_coordinate_2,\
                            *args):
    ###################################################
    ## Calcualtes the energy of a individual particles
    ## from coordinates of it edges.
    ## Note that edge numbered with "0" must be red.
    ###################################################
    ###############
    (k_red, k_blue, k_soft, \
     k_area_red, k_area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     Area0_red, Area0_blue, \
     epsilon, width)= args
    ###############
    energy_sum_temp = 0.0
    ## spring energy
    energy_sum_spring_temp = 0.0
    for ind_sp in range(3):
        ## Red triangles
        length_x = edge_coordinate_1[ind_sp, 0] - edge_coordinate_1[(ind_sp+1)%3, 0]
        length_y = edge_coordinate_1[ind_sp, 1] - edge_coordinate_1[(ind_sp+1)%3, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = 0.5*k_red*(lengthSQ - l0_red)*(lengthSQ - l0_red)
        energy_sum_spring_temp = energy_sum_spring_temp + en_temptemp
        ## Blue triangles
        length_x = edge_coordinate_2[ind_sp, 0] - edge_coordinate_2[(ind_sp+1)%3, 0]
        length_y = edge_coordinate_2[ind_sp, 1] - edge_coordinate_2[(ind_sp+1)%3, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = 0.5*k_blue*(lengthSQ - l0_blue)*(lengthSQ - l0_blue)
        energy_sum_spring_temp = energy_sum_spring_temp + en_temptemp
        ## soft springs
        length_x = edge_coordinate_1[ind_sp, 0] - edge_coordinate_2[(ind_sp+1)%3, 0]
        length_y = edge_coordinate_1[ind_sp, 1] - edge_coordinate_2[(ind_sp+1)%3, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = 0.5*k_soft*(lengthSQ - l0_soft)*(lengthSQ - l0_soft)
        length_x = edge_coordinate_1[ind_sp, 0] - edge_coordinate_2[(ind_sp+2)%3, 0]
        length_y = edge_coordinate_1[ind_sp, 1] - edge_coordinate_2[(ind_sp+2)%3, 1]
        lengthSQ = math.sqrt(length_x*length_x + length_y*length_y)
        en_temptemp = en_temptemp + \
                      0.5*k_soft*(lengthSQ - l0_soft)*(lengthSQ - l0_soft)
        energy_sum_spring_temp = energy_sum_spring_temp + en_temptemp
    ## Area term in energy
    area_red = FindAreaOfTriangles(edge_coordinate_1)
    area_blue = FindAreaOfTriangles(edge_coordinate_2)
    if area_red<0 or area_blue<0:
        area_red = 10000000.0
        area_blue = 10000000.0
    energy_sum_temp = energy_sum_spring_temp + \
                      0.5*(k_area_red/Area0_red)*\
                      (area_red - Area0_red)*\
                      (area_red - Area0_red)+\
                       + 0.5*(k_area_blue/Area0_blue)*\
                      (area_blue - Area0_blue)*\
                      (area_blue - Area0_blue)
    #########
    return energy_sum_temp
########################################################
########################################################

########################################################
########################################################
def FunctionBulkEnergy(l_bulk_temp, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_area_red, k_area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     Area0_red, Area0_blue, \
     epsilon)= args
    ###############
    Area_bulk_temp = (math.sqrt(3.0)/4.0)*\
                     l_bulk_temp*l_bulk_temp
    energy_temp = (3.0/2.0)*k_red*\
                  (l_bulk_temp - l0_red)*\
                  (l_bulk_temp - l0_red) + \
                  (3.0/2.0)*k_blue*\
                  (l_bulk_temp - l0_blue)*\
                  (l_bulk_temp - l0_blue) + \
                  (1.0/2.0)*(k_area_red/Area0_red) *\
                  (Area_bulk_temp - Area0_red)*\
                  (Area_bulk_temp - Area0_red) + \
                  (1.0/2.0)*(k_area_blue/Area0_blue) *\
                  (Area_bulk_temp - Area0_blue)*\
                  (Area_bulk_temp - Area0_blue)+ \
                  (6.0/2.0)*k_soft*\
                  (l_bulk_temp - l0_soft)*\
                  (l_bulk_temp - l0_soft)
    energy_temp = energy_temp
    ##
    return energy_temp
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
    res = fmin_cg(FunctionBulkEnergy, l0_0, \
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
