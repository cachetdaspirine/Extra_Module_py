########################################################
## Functions which are specific to the model ###########
########################################################
import time
import numpy as np
from scipy.optimize import fmin_cg
from scipy.optimize import minimize
import math
########################################################
########################################################

########################################################
########################################################
from ModelFunctions_v2 import *
########################################################
########################################################

########################################################
########################################################
def XYfromDOF_SecondConfiguation(xy_dof):
    ## Determine xy values from dof
    width = int((np.size(xy_dof, 0) - 9)/8)
    N_edges = 3*(2*width + 2)
    ##
    xy_list = np.zeros([N_edges, 2])
    xy_list[2, 0] = xy_dof[np.size(xy_dof, 0)-1]
    deltaX_temp = xy_dof[np.size(xy_dof, 0)-1] - xy_dof[0]
    for ind_d in range(N_edges):
        ind_dof = ind_d - (ind_d//3)
        if ind_d%3==2:
            if ind_d>2:
                x_temp = deltaX_temp + xy_list[ind_d-2, 0]
                xy_list[ind_d, 0] = x_temp
                xy_list[ind_d, 1] = xy_list[ind_d-2, 1]
        else:
            xy_list[ind_d, 0] = xy_dof[2*ind_dof]
            xy_list[ind_d, 1] = xy_dof[2*ind_dof+1]
    ##
    ##
    return xy_list
########################################################
########################################################

########################################################
########################################################
def DOFfromXY_SecondConfiguation(xy_list):
    ## Determine DOF values from edge coordinates
    N_edges = np.size(xy_list, 0)
    width = int((N_edges-6)/6)
    N_dof = int(8*width + 9)
    ##
    xy_dof = np.zeros(N_dof)
    for ind_d in range(N_edges):
        ind_dof = ind_d - (ind_d//3)
        if ind_d%3!=2:
            xy_dof[2*ind_dof] = xy_list[ind_d, 0]
            xy_dof[2*ind_dof+1] = xy_list[ind_d, 1]
    ##
    xy_dof[np.size(xy_dof, 0)-1] = xy_list[2, 0]
    ##
    return xy_dof
########################################################
########################################################

########################################################
########################################################
def DetermineListOfEdgesOfHexagons_SecondConfiguation(width):
    ## List of edges belongs to that hexagon
    N_hex = 2*width
    hexagons_edges = np.zeros([N_hex, 6],dtype=int)
    for ind_h in range(N_hex):
        if ind_h%2==0:
            ind_0 = (ind_h/2)*6
        elif ind_h%2==1:
            ind_0 = ((ind_h-1)/2)*6 + 4
        hexagons_edges[ind_h, 0] = ind_0
        hexagons_edges[ind_h, 1] = ind_0 + 1
        hexagons_edges[ind_h, 2] = ind_0 + 4
        hexagons_edges[ind_h, 3] = ind_0 + 7
        hexagons_edges[ind_h, 4] = ind_0 + 6
        hexagons_edges[ind_h, 5] = ind_0 + 3
    return hexagons_edges

########################################################
########################################################
def compute_energy_along_width_type_2(xy_dof,*args):
    Mc,rho0 = args
    xy_per_hexagon = XYfromDOF_SecondConfiguation(xy_dof)
    width = (xy_per_hexagon.shape[0] - 6)//6
    hexagons_edges_index = DetermineListOfEdgesOfHexagons_SecondConfiguation(width)
    E = np.zeros(2*width,dtype=float)
    for w in range(2*width):
        edge_coordinate = np.zeros([6, 2])
        #edge_coordinate2 = np.zeros([6, 2])
        for i in range(6):
            ind_edge = hexagons_edges_index[w,i]
            edge_coordinate[i, 0] = xy_per_hexagon[ind_edge, 0]
            edge_coordinate[i, 1] = xy_per_hexagon[ind_edge, 1]

            #ind_edge = hexagons_edges_index[2*w+1,i]
            #edge_coordinate2[i, 0] = xy_per_hexagon[ind_edge, 0]
            #edge_coordinate2[i, 1] = xy_per_hexagon[ind_edge, 1]
        #E[w] = np.mean([CalculateParticleEnergy(edge_coordinate,*args),CalculateParticleEnergy(edge_coordinate2,*args)])
        E[w] = CalculateParticleEnergy(edge_coordinate,*args)
    return E
########################################################
########################################################
def CalcualteTotalEnergy_SecondConfiguation(xy_dof, *args):
    ## Calculates total energy from the values
    ## of degrees of freedom
    ################################################
    ###############
    (M_ij, rho0_vec)= args
    ################################################
    ## Frist determine the edge coordinates from DOF
    xy_list = XYfromDOF_SecondConfiguation(xy_dof)
    ################################################
    width = int((np.size(xy_list, 0) - 6)/6)
    N_hex = 2*width
    ################################################
    ## List of edges belongs to that hexagon
    hexagons_edges = DetermineListOfEdgesOfHexagons_SecondConfiguation(width)
    ################################################
    edge_coordinate = np.zeros([6, 2])
    energyTotal_temp = 0.0
    for ind_hg_i in range(N_hex):
        for ind_eg in range(6):
            ind_edge_0 = hexagons_edges[ind_hg_i, ind_eg]
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
def GradDOFfromXY_SecondConfiguation(df_array_temp, df_array_temp_1d):
    N_dof = np.size(df_array_temp_1d, 0)
    N_edges = np.size(df_array_temp, 0)
    for ind_df in range(N_edges):
        if ind_df==2:
            df_array_temp_1d[N_dof-1] = df_array_temp_1d[N_dof-1] + \
                                        df_array_temp[ind_df, 0]
            df_array_temp_1d[1] = df_array_temp_1d[1] + \
                                  df_array_temp[ind_df, 1]
        elif ind_df%3!=2:
            ind_q = (2*(ind_df//3) + ind_df%3)*2
            df_array_temp_1d[ind_q] = df_array_temp_1d[ind_q] + \
                                      df_array_temp[ind_df, 0]
            df_array_temp_1d[ind_q+1] = df_array_temp_1d[ind_q+1] + \
                                        df_array_temp[ind_df, 1]
        elif ind_df%3==2 and ind_df!=2:
            df_array_temp_1d[0] = df_array_temp_1d[0] - \
                                  df_array_temp[ind_df, 0]
            df_array_temp_1d[N_dof-1] = df_array_temp_1d[N_dof-1] \
                                        + df_array_temp[ind_df, 0]
            ind_q = (4*(ind_df//3))
            df_array_temp_1d[ind_q] = df_array_temp_1d[ind_q] + \
                                      df_array_temp[ind_df, 0]
            df_array_temp_1d[ind_q+1] = df_array_temp_1d[ind_q+1] + \
                                        df_array_temp[ind_df, 1]
    return df_array_temp_1d
########################################################
########################################################

########################################################
########################################################
def gradf_f_SecondConfiguation(xy_dof, *args):
    ###############
    (M_ij, rho0_vec)= args
    ###############
    ##
    nodes_coordinates = XYfromDOF_SecondConfiguation(xy_dof)
    width_temp = int((np.size(nodes_coordinates, 0) - 6)/6)
    ##
    hexagons_nodes = DetermineListOfEdgesOfHexagons_SecondConfiguation(width_temp)
    ##
    N_nodes = np.size(nodes_coordinates, 0)
    N_hexagons = np.size(hexagons_nodes, 0)
    df_array_temp = np.zeros([N_nodes, 2])
    edge_coordinate = np.zeros([6, 2])
    for ind_df in range(N_hexagons):
        for ind_eg in range(6):
            ind_edge_0 = hexagons_nodes[ind_df, ind_eg]
            edge_coordinate[ind_eg, 0] = nodes_coordinates[ind_edge_0, 0]
            edge_coordinate[ind_eg, 1] = nodes_coordinates[ind_edge_0, 1]
        force_unit_temp = ForcesUnit(edge_coordinate, *args)
        for ind_df_j in range(6):
            ind_df_temp = hexagons_nodes[ind_df, ind_df_j]
            df_array_temp[ind_df_temp, 0] = df_array_temp[ind_df_temp, 0]\
                                            + force_unit_temp[ind_df_j, 0]
            df_array_temp[ind_df_temp, 1] = df_array_temp[ind_df_temp, 1]\
                                            + force_unit_temp[ind_df_j, 1]
    ##
    df_array_temp_1d = np.zeros(np.size(xy_dof, 0))
    df_array_temp_1d = GradDOFfromXY_SecondConfiguation(df_array_temp, df_array_temp_1d)
    ## Normalize to per particle values
    for ind_n in range(np.size(df_array_temp_1d, 0)):
        df_array_temp_1d[ind_n] = df_array_temp_1d[ind_n]/N_hexagons
    ##
    return (df_array_temp_1d)
########################################################
########################################################

########################################################
########################################################
def InitialXYlist_SecondConfiguation(width):
    N_edges = 3*(2*width + 2)
    l0_0 = 1.0
    l0_soft = l0_0*math.sqrt((1.0/3.0))
    ##
    xy_list = np.zeros([N_edges,2])
    for ind_d in range(N_edges):
        ind_x = ind_d//3
        ind_y = ind_d%3
        y_temp = ind_x*(l0_0/2.0)
        if ind_x%2==0:
            x_ini_list = [1.0/2.0, 3.0/2.0, 7.0/2.0]
            x_temp = l0_soft*x_ini_list[ind_y]
        elif ind_x%2==1:
            x_ini_list = [0.0, 2.0, 3.0]
            x_temp = l0_soft*x_ini_list[ind_y]
        ##
        xy_list[ind_d, 0] = x_temp
        xy_list[ind_d, 1] = y_temp
        ##
    return xy_list
########################################################
########################################################

########################################################
########################################################
def FunctionFilamentEnergy_SecondConfiguation(width, *args):
    ###############
    (M_ij, rho0_vec)= args
    ###############
    xy_list_0 = InitialXYlist_SecondConfiguation(width)
    xy_dof_0 = DOFfromXY_SecondConfiguation(xy_list_0)
    ###############
    #res = fmin_cg(CalcualteTotalEnergy_SecondConfiguation, xy_dof_0, \
    #              fprime=gradf_f_SecondConfiguation, \
    #              args=args, full_output=1, disp=0)
    res = minimize(CalcualteTotalEnergy_SecondConfiguation,xy_dof_0, \
                    jac = gradf_f_SecondConfiguation, \
                    args=args,method='BFGS')
    res_min = res.x
    ###############
    en_w_temp = CalcualteTotalEnergy_SecondConfiguation(res_min, *args)
    ###############
    return en_w_temp
########################################################
########################################################
def get_E_along_width_type_2(width,*args):
    (M_ij, rho0_vec)= args
    ###############
    xy_list_0 = InitialXYlist_SecondConfiguation(width)
    xy_dof_0 = DOFfromXY_SecondConfiguation(xy_list_0)
    ###############
    #res = fmin_cg(CalcualteTotalEnergy_FirstConfiguation, xy_dof_0, \
    #              fprime=gradf_f_FirstConfiguation, \
    #              args=args, gtol=10.0**(-5), full_output=1, disp=0)
    res = minimize(CalcualteTotalEnergy_SecondConfiguation,xy_dof_0, \
                    jac = gradf_f_SecondConfiguation, \
                    args=args,method='BFGS')
    if not res.success:
        print(res.message)
    #res_min = res[0]
    #if zero_curvature == True:
    #    res.x[-1] = 0
    E = compute_energy_along_width_type_2(res.x, *args)
    #if full==True:
    #    return E
    if E.shape[0]%2 == 0:
        return np.mean([E[:E.shape[0]//2],np.flip(E[E.shape[0]//2:])],axis=0)
    else:
        return np.append(np.mean([E[:E.shape[0]//2],np.flip(E[E.shape[0]//2+1:])],axis=0),E[E.shape[0]//2])
########################################################
########################################################
