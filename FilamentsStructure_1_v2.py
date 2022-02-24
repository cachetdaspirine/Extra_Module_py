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
def compute_energy_along_width(xy_dof,*args):
    Mc,rho0 = args
    xy_per_hexagon = XYfromDOF_FirstConfiguation(xy_dof)
    width = (xy_per_hexagon.shape[0] - 4)//4
    hexagons_edges_index = DetermineListOfEdgesOfHexagons_FirstConfiguation(width)
    E = np.zeros(width,dtype=float)
    for w in range(width):
        edge_coordinate = np.zeros([6, 2])
        for i in range(6):
            ind_edge = hexagons_edges_index[w,i]
            edge_coordinate[i, 0] = xy_per_hexagon[ind_edge, 0]
            edge_coordinate[i, 1] = xy_per_hexagon[ind_edge, 1]
        E[w] = CalculateParticleEnergy(edge_coordinate,*args)
    return E

def CalcualteTotalEnergy_FirstConfiguation(xy_dof, *args):
    ## Calculates total energy from the values
    ## of degrees of freedom
    ################################################
    ###############
    (M_ij, rho0_vec)= args
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
    (M_ij, rho0_vec)= args
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
    df_array_temp_1d = GradDOFfromXY_FirstConfiguation(df_array_temp, df_array_temp_1d, xy_dof)
    ## Normalize to per particle values
    for ind_n in range(np.size(df_array_temp_1d, 0)):
        df_array_temp_1d[ind_n] = df_array_temp_1d[ind_n]/N_hexagons
    ##
##    df_n = NumericalGradient(xy_dof, *args)
##    fark_st = np.linalg.norm(df_array_temp_1d - df_n)
##    if fark_st>0.01:
##        print('df diff: ', fark_st)
##        for ind_pp in range(np.size(xy_dof, 0)):
##            print('%.10f\n%.10f\n%.10f\n' %(df_array_temp_1d[ind_pp] - df_n[ind_pp], \
##                                            df_array_temp_1d[ind_pp], df_n[ind_pp]))
##
    ##
    return (df_array_temp_1d)
########################################################
########################################################

########################################################
########################################################
def DetermineListOfEdgesOfHexagons_FirstConfiguation(width):
    ## List of edges belongs to that hexagon
    N_hex = width
    hexagons_edges = np.zeros([N_hex, 6],dtype=int)
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
def InitialXYlist_FirstConfiguation(width):
    N_edges = 4*width + 4
    curvature_temp = 10.0**(-4) ## radius of curvature
    #curvature_temp = ((1-eps)/2*np.sqrt(1-3*eps**2)/(np.sqrt(3)*eps))**(-1)
    l0_0 = 1.0
    l0_soft = l0_0*math.sqrt((1.0/3.0))
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
    (M_ij, rho0_vec)= args
    ###############
    (xy_list_0, C_ini) = InitialXYlist_FirstConfiguation(width)
    xy_dof_0 = DOFfromXY_FirstConfiguation(xy_list_0, C_ini)
    ###############
    #res = fmin_cg(CalcualteTotalEnergy_FirstConfiguation, xy_dof_0, \
    #              fprime=gradf_f_FirstConfiguation, \
    #              args=args, gtol=10.0**(-5), full_output=1, disp=0)
    res = minimize(CalcualteTotalEnergy_FirstConfiguation,xy_dof_0, \
                    jac = gradf_f_FirstConfiguation, \
                    args=args,method='BFGS')
    #res_min = res[0]
    #res_min = res.x
    ###############
    #df_0 = gradf_f_FirstConfiguation(xy_dof_0, *args)
    #df_0_n = NumericalGradient(xy_dof_0, *args)
##    print('initial grad diff: ', np.linalg.norm(df_0 - df_0_n), \
##          sum(abs(df_0)))
    #xy_fin = res[0]
    xy_fin = res.x
#    df_r = gradf_f_FirstConfiguation(xy_fin, *args)
##    print('final grad total: ', sum(abs(df_r)))
#    Rcuv = xy_fin[np.size(xy_fin,0)-1]
    if not res.success:
        print(res.message)
    en_w_temp = CalcualteTotalEnergy_FirstConfiguation(xy_fin, *args)
    ###############
    return en_w_temp
########################################################
########################################################
def get_E_along_width(width,*args):
    (M_ij, rho0_vec)= args
    ###############
    (xy_list_0, C_ini) = InitialXYlist_FirstConfiguation(width)
    xy_dof_0 = DOFfromXY_FirstConfiguation(xy_list_0, C_ini)
    ###############
    #res = fmin_cg(CalcualteTotalEnergy_FirstConfiguation, xy_dof_0, \
    #              fprime=gradf_f_FirstConfiguation, \
    #              args=args, gtol=10.0**(-5), full_output=1, disp=0)
    res = minimize(CalcualteTotalEnergy_FirstConfiguation,xy_dof_0, \
                    jac = gradf_f_FirstConfiguation, \
                    args=args,method='BFGS')
    if not res.success:
        print(res.message)
    #res_min = res[0]
    E = compute_energy_along_width(res.x, *args)
    if E.shape[0]%2 == 0:
        return np.mean([E[:E.shape[0]//2],np.flip(E[E.shape[0]//2:])],axis=0)
    else:
        return np.append(np.mean([E[:E.shape[0]//2],np.flip(E[E.shape[0]//2+1:])],axis=0),E[E.shape[0]//2])

########################################################
########################################################
def NumericalGradient(xy_dof, *args):
    ###############
    (M_ij, rho0_vec)= args
    ###############
    x_delta = 10.0**(-6)
    ###############
    energy_1 = CalcualteTotalEnergy_FirstConfiguation(xy_dof, *args)
    df_1 = np.zeros(np.size(xy_dof, 0))
    for ind_n in range(np.size(xy_dof, 0)):
        xy_temp = np.zeros(np.size(xy_dof, 0))
        for ind_nj in range(np.size(xy_dof, 0)):
            xy_temp[ind_nj] = xy_dof[ind_nj]
        xy_temp[ind_n] = xy_temp[ind_n] + x_delta
        energy_2 = CalcualteTotalEnergy_FirstConfiguation(xy_temp, *args)
        df_1[ind_n] = (energy_2 - energy_1)/x_delta
    return df_1
########################################################
########################################################


########################################################
########################################################
