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
from ModelFunctionsHex_v1 import *
########################################################
########################################################

########################################################
########################################################
def CalcualteTotalEnergy_FirstConfiguation(xy_dof, *args):
    ## Calculates total energy from the values
    ## of degrees of freedom
    ################################################
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ################################################
    ## Frist determine the edge coordinates from DOF
    xy_list = XYfromDOF_FirstConfiguation(xy_dof)
    ################################################
    width = int((np.size(xy_list, 0) - 4)/4)
    N_hex = width
    ################################################
    ## List of edges belongs to that hexagon
    hexagons_edges = DetermineListOfEdgesOfHexagons_FirstConfiguation(width)
    ################################################
    edge_coordinate = np.zeros([6, 2])
    energyTotal_temp = 0.0
    for ind_hg_i in range(N_hex):
        for ind_eg in range(6):
            ind_edge_0 = int(hexagons_edges[ind_hg_i, ind_eg])
            edge_coordinate[ind_eg, 0] = xy_list[ind_edge_0, 0]
            edge_coordinate[ind_eg, 1] = xy_list[ind_edge_0, 1]
        ##
        energy_unit_sum_temp = CalculateParticleEnergy(edge_coordinate, \
                                                       *args)
        ##
        energyTotal_temp = energyTotal_temp + energy_unit_sum_temp
    ##
    energyTotal = energyTotal_temp/N_hex
    #################################################
    return energyTotal
########################################################
########################################################

########################################################
########################################################
def GradDOFfromXY_FirstConfiguation(df_array_temp, df_array_temp_1d, xy_dof):
    N_dof = np.size(df_array_temp_1d, 0)
    N_edges = np.size(df_array_temp, 0)
    ##
    curvature_temp = xy_dof[N_dof-1]
    ##
    xy_dof_temp = np.zeros(N_dof)
    for ind_df in range(N_dof-1):
        xy_dof_temp[ind_df] = xy_dof[ind_df]
    xy_dof_temp[N_dof-1] = 0.0
    ##
    xy0_list = XYfromDOF_FirstConfiguation(xy_dof_temp)
    ##
    for ind_df in range(N_edges):
        x_temp = xy0_list[ind_df, 0]
        y_temp = xy0_list[ind_df, 1]
        ####
        phi_temp = curvature_temp*x_temp
        r_temp = 1.0/curvature_temp
        ####
        SinPhi = math.sin(phi_temp)
        CosPhi = math.cos(phi_temp)
        ####
        coef1 = CosPhi*(r_temp + y_temp)/r_temp
        coef2 = SinPhi
        coef3 = - SinPhi*(r_temp + y_temp)/r_temp
        coef4 = CosPhi
        ####
        if ind_df%2==0:
            ##
            df_array_temp_1d[ind_df] = df_array_temp_1d[ind_df] + \
                                       coef1*df_array_temp[ind_df, 0]
            df_array_temp_1d[ind_df+1] = df_array_temp_1d[ind_df+1] + \
                                         coef2*df_array_temp[ind_df, 0]
            ##
            df_array_temp_1d[ind_df] = df_array_temp_1d[ind_df] + \
                                       coef3*df_array_temp[ind_df, 1]
            df_array_temp_1d[ind_df+1] = df_array_temp_1d[ind_df+1] + \
                                         coef4*df_array_temp[ind_df, 1]
        elif ind_df==1:
            ##
            df_array_temp_1d[N_dof-3] = df_array_temp_1d[N_dof-3] + \
                                        coef1*df_array_temp[ind_df, 0]
            df_array_temp_1d[N_dof-2] = df_array_temp_1d[N_dof-2] + \
                                        coef2*df_array_temp[ind_df, 0]
            ##
            df_array_temp_1d[N_dof-3] = df_array_temp_1d[N_dof-3] + \
                                        coef3*df_array_temp[ind_df, 1]
            df_array_temp_1d[N_dof-2] = df_array_temp_1d[N_dof-2] + \
                                        coef4*df_array_temp[ind_df, 1]
        elif ind_df%2==1 and ind_df>1:
            ##
            df_array_temp_1d[0] = df_array_temp_1d[0] - \
                                  coef1*df_array_temp[ind_df, 0]
            df_array_temp_1d[N_dof-3] = df_array_temp_1d[N_dof-3] + \
                                        coef1*df_array_temp[ind_df, 0]
            df_array_temp_1d[ind_df-1] = df_array_temp_1d[ind_df-1] + \
                                         coef1*df_array_temp[ind_df, 0]
            df_array_temp_1d[1] = df_array_temp_1d[1] - \
                                  coef2*df_array_temp[ind_df, 0]
            df_array_temp_1d[N_dof-2] = df_array_temp_1d[N_dof-2] + \
                                        coef2*df_array_temp[ind_df, 0]
            df_array_temp_1d[ind_df] = df_array_temp_1d[ind_df] + \
                                       coef2*df_array_temp[ind_df, 0]
            ##
            df_array_temp_1d[0] = df_array_temp_1d[0] - \
                                  coef3*df_array_temp[ind_df, 1]
            df_array_temp_1d[N_dof-3] = df_array_temp_1d[N_dof-3] + \
                                        coef3*df_array_temp[ind_df, 1]
            df_array_temp_1d[ind_df-1] = df_array_temp_1d[ind_df-1] + \
                                         coef3*df_array_temp[ind_df, 1]
            df_array_temp_1d[1] = df_array_temp_1d[1] - \
                                  coef4*df_array_temp[ind_df, 1]
            df_array_temp_1d[N_dof-2] = df_array_temp_1d[N_dof-2] + \
                                        coef4*df_array_temp[ind_df, 1]
            df_array_temp_1d[ind_df] = df_array_temp_1d[ind_df] + \
                                       coef4*df_array_temp[ind_df, 1]
        ## Derivtive Terms of Curvature
        coef_C1 = - r_temp*r_temp*SinPhi \
                  + (r_temp + y_temp)*x_temp*CosPhi
        coef_C2 = - r_temp*r_temp*CosPhi + r_temp*r_temp \
                  - (r_temp + y_temp)*x_temp*SinPhi
        df_array_temp_1d[N_dof-1] = df_array_temp_1d[N_dof-1] + \
                                    coef_C1*df_array_temp[ind_df, 0]
        df_array_temp_1d[N_dof-1] = df_array_temp_1d[N_dof-1] + \
                                    coef_C2*df_array_temp[ind_df, 1]
    return df_array_temp_1d
########################################################
########################################################

########################################################
########################################################
def gradf_f_FirstConfiguation(xy_dof, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    ##
    nodes_coordinates = XYfromDOF_FirstConfiguation(xy_dof)
    width_temp = int((np.size(nodes_coordinates, 0) - 4)/4)
    ##
    hexagons_nodes = DetermineListOfEdgesOfHexagons_FirstConfiguation(width_temp)
    ##
    N_nodes = np.size(nodes_coordinates, 0)
    N_hexagons = np.size(hexagons_nodes, 0)
    df_array_temp = np.zeros([N_nodes, 2])
    edge_coordinate = np.zeros([6, 2])
    for ind_df in range(N_hexagons):
        for ind_eg in range(6):
            ind_edge_0 = int(hexagons_nodes[ind_df, ind_eg])
            edge_coordinate[ind_eg, 0] = nodes_coordinates[ind_edge_0, 0]
            edge_coordinate[ind_eg, 1] = nodes_coordinates[ind_edge_0, 1]
        force_unit_temp = ForcesUnit(edge_coordinate, *args)
        for ind_df_j in range(6):
            ind_df_temp = int(hexagons_nodes[ind_df, ind_df_j])
            df_array_temp[ind_df_temp, 0] = df_array_temp[ind_df_temp, 0]\
                                            + force_unit_temp[ind_df_j, 0]
            df_array_temp[ind_df_temp, 1] = df_array_temp[ind_df_temp, 1]\
                                            + force_unit_temp[ind_df_j, 1]
    ##
    df_array_temp_1d = np.zeros(np.size(xy_dof, 0))
    df_array_temp_1d = GradDOFfromXY_FirstConfiguation(df_array_temp, df_array_temp_1d, xy_dof)
    ## Normalize to per particle values
    for ind_n in range(np.size(df_array_temp_1d, 0)):
        df_array_temp_1d[ind_n] = df_array_temp_1d[ind_n]/N_hexagons
    ##
    return (df_array_temp_1d)
########################################################
########################################################

########################################################
########################################################
def DetermineListOfEdgesOfHexagons_FirstConfiguation(width):
    ## List of edges belongs to that hexagon
    N_hex = width
    hexagons_edges = np.zeros([N_hex, 6])
    for ind_h in range(N_hex):
        if ind_h%2==0:
            ind_0 = (ind_h)*4
            hexagons_edges[ind_h, 0] = ind_0
            hexagons_edges[ind_h, 1] = ind_0 + 3
            hexagons_edges[ind_h, 2] = ind_0 + 5
            hexagons_edges[ind_h, 3] = ind_0 + 6
            hexagons_edges[ind_h, 4] = ind_0 + 4
            hexagons_edges[ind_h, 5] = ind_0 + 2
        elif ind_h%2==1:
            ind_0 = (ind_h-1)*4 + 5
            hexagons_edges[ind_h, 0] = ind_0
            hexagons_edges[ind_h, 1] = ind_0 + 2
            hexagons_edges[ind_h, 2] = ind_0 + 4
            hexagons_edges[ind_h, 3] = ind_0 + 6
            hexagons_edges[ind_h, 4] = ind_0 + 3
            hexagons_edges[ind_h, 5] = ind_0 + 1
    return hexagons_edges
########################################################
########################################################

########################################################
########################################################
def BendTheFilament_FirstConfiguation(xy_nodes_array_temp, curvature_temp):
    N_edges = np.size(xy_nodes_array_temp, 0)
    ## Bending of the filament.
    for ind_d in range(N_edges):
        x_temp = xy_nodes_array_temp[ind_d, 0]
        y_temp = xy_nodes_array_temp[ind_d, 1]
        ####
        phi_temp = curvature_temp * x_temp
        r_temp = 1.0/curvature_temp
        ####
        SinPhi = math.sin(phi_temp)
        CosPhi = math.cos(phi_temp)
        ####
        xy_nodes_array_temp[ind_d, 0] = (r_temp + y_temp)*\
                                        SinPhi
        xy_nodes_array_temp[ind_d, 1] = (r_temp + y_temp)*\
                                        CosPhi \
                                        - r_temp
    return xy_nodes_array_temp
########################################################
########################################################

########################################################
########################################################
def DOFfromXY_FirstConfiguation(xy_list, curvature_temp):
    ## Determine DOF values from edge coordinates
    N_edges = np.size(xy_list, 0)
    width = int((N_edges-4)/4)
    N_dof = int(4*width + 4 + 3)
    ##
    xy_dof = np.zeros(N_dof)
    for ind_d in range(N_edges):
        ind_dof = (ind_d//2)
        if ind_d%2==0:
            xy_dof[2*ind_dof] = xy_list[ind_d, 0]
            xy_dof[2*ind_dof+1] = xy_list[ind_d, 1]
    ##
    xy_dof[np.size(xy_dof, 0)-3] = xy_list[1, 0]
    xy_dof[np.size(xy_dof, 0)-2] = xy_list[1, 1]
    xy_dof[np.size(xy_dof, 0)-1] = curvature_temp
    ##
    return xy_dof
########################################################
########################################################

########################################################
########################################################
def XYfromDOF_FirstConfiguation(xy_dof):
    ## Determine xy values from dof
    N_dof = np.size(xy_dof, 0)
    width = int((N_dof-7)/4)
    N_edges = int(4*width + 4)
    ##
    curvature_temp = xy_dof[N_dof-1]
    ##
    xy_list = np.zeros([N_edges, 2])
    ##
    xy_list[1, 0] = xy_dof[N_dof-3]
    xy_list[1, 1] = xy_dof[N_dof-2]
    ##
    delta_x_temp = xy_dof[N_dof-3] - xy_dof[0]
    delta_y_temp = xy_dof[N_dof-2] - xy_dof[1]
    for ind_d in range(N_edges):
        ind_dof = (ind_d//2)
        if ind_d%2==0:
            ind_dof = ind_d
            xy_list[ind_d, 0] = xy_dof[ind_dof]
            xy_list[ind_d, 1] = xy_dof[ind_dof + 1]
        elif ind_d%2==1 and ind_d>1:
            ind_dof = ind_d-1
            xy_list[ind_d, 0] = delta_x_temp + xy_dof[ind_dof]
            xy_list[ind_d, 1] = delta_y_temp + xy_dof[ind_dof + 1]
    ##
    if curvature_temp!=0.0:
        xy_list = BendTheFilament_FirstConfiguation(xy_list, curvature_temp)
    ##
    return xy_list
########################################################
########################################################

########################################################
########################################################
def InitialXYlist_FirstConfiguation(width, l0_0, l0_soft):
    N_edges = 4*width + 4
    curvature_temp = 0.001 ## radius of curvature
    ##
    xy_list = np.zeros([N_edges,2])
    for ind_d in range(N_edges):
        ind_l = ind_d//2
        ind_y = ind_d%2
        if ind_l%4==0 or ind_l%4==3:
            xy_list[ind_d, 0] = (ind_y + 0.5)*l0_0
        else:
            xy_list[ind_d, 0] = (ind_y)*l0_0
        ##
        if ind_l%2==0:
            xy_list[ind_d, 1] = (ind_l/2.0)*(3.0/2.0)*l0_soft
        elif ind_l%2==1:
            xy_list[ind_d, 1] = (((ind_l-1.0)/2.0)*(3.0/2.0)+1.0/2.0)*l0_soft
        ##
    #xy_list = BendTheFilament(xy_list, curvature_temp)
    ##
    return (xy_list, curvature_temp)
########################################################
########################################################

########################################################
########################################################
def FunctionFilamentEnergy_FirstConfiguation(width, *args):
    ###############
    (k_red, k_blue, k_soft, \
     k_Area_red, k_Area_blue, \
     l0_0, l0_red, l0_blue, l0_soft, \
     area_0_red, area_0_blue, \
     epsilon)= args
    ###############
    (xy_list_0, C_ini) = InitialXYlist_FirstConfiguation(width, l0_0, l0_soft)
    xy_dof_0 = DOFfromXY_FirstConfiguation(xy_list_0, C_ini)
    ###############
    res = fmin_cg(CalcualteTotalEnergy_FirstConfiguation, xy_dof_0, \
                  fprime=gradf_f_FirstConfiguation, \
                  args=args, full_output=1, disp=0)
    res_min = res[0]
    ###############
    if res[4]!=0:
        print('\n')
        print('Error: There is an error in the filament energy calculations in the first sturcture')
        print('\n')
        print(res[4])
    ###############
    en_w_temp = CalcualteTotalEnergy_FirstConfiguation(res_min, *args)
    ###############
    return en_w_temp
########################################################
########################################################

########################################################
########################################################
