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
#from mapping import *
mapping = np.array([[0, 2, 4, 6, 8,10, 0, 4, 8],\
                    [4, 6, 8,10, 0, 2, 2, 6,10]])
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
    (M_ij, rho0_vec)= args
    ###############
    ##
    q_vec = EdgeCoordsToQVEC(edge_coordinate)
    rho_vec = DetermineLocalCoordinates(q_vec)
    ##
    ## Note that rho_vec is the square of the lengths
    ## and rho0_vec is rest lengths, not squared.
    ##
    # eta_vec = np.zeros(9)
    # for ind_i in range(9):
    #     eta_vec[ind_i] = (rho_vec[ind_i]/rho0_vec[ind_i] \
    #                       - rho0_vec[ind_i])\
    #                       /(2.0)
    eta_vec = (rho_vec/rho0_vec -rho0_vec)*0.5
    ##
    # energy_t = 0.0
    # for ind_i in range(9):
    #     for ind_j in range(9):
    #         energy_t = energy_t \
    #                    + (1.0/2.0)*eta_vec[ind_i]*\
    #                    M_ij[ind_i, ind_j]*\
    #                    eta_vec[ind_j]
    energy_t = 0.5*np.dot(eta_vec,np.dot(M_ij,eta_vec))
    ##
    #########
    return energy_t
########################################################
########################################################

########################################################
########################################################
def EdgeCoordsToQVEC(edge_coordinate):
    q_vec = np.zeros(12)
    for ind_q in range(6):
        q_vec[2*ind_q] = edge_coordinate[ind_q, 0]
        q_vec[2*ind_q+1] = edge_coordinate[ind_q, 1]
    return q_vec
########################################################
########################################################

########################################################
########################################################
def DetermineLocalCoordinates(q_vec):
    ##
    #mapping = UseMapping()
    ##
    rho_vec = np.zeros(9)
    for ind_i in range(9):
        ind_1 = mapping[0, ind_i]
        ind_2 = mapping[1, ind_i]
        rho_vec[ind_i] = (q_vec[ind_1] - q_vec[ind_2])\
                         *(q_vec[ind_1] - q_vec[ind_2])\
                         + (q_vec[ind_1+1] - q_vec[ind_2+1])\
                         *(q_vec[ind_1+1] - q_vec[ind_2+1])
    return rho_vec
########################################################
########################################################

########################################################
########################################################
def ForcesUnit(edge_coordinate, *args):
    ###############
    (M_ij, rho0_vec)= args
    ###############
    force_unit_temp = np.zeros([6, 2])
    ##
    #mapping = UseMapping()
    ##
    q_vec = EdgeCoordsToQVEC(edge_coordinate)
    ##
    DqDeta_ij = np.zeros([9,12])
    for ind_i in range(9):
        ind1 = mapping[0, ind_i]
        ind2 = mapping[1, ind_i]
        ##
        DqDeta_ij[ind_i, ind1] += (q_vec[ind1]-q_vec[ind2])\
                                    /rho0_vec[ind_i]

        ##
        DqDeta_ij[ind_i, ind2] += -(q_vec[ind1]-q_vec[ind2])\
                                    /rho0_vec[ind_i]
        ##
        DqDeta_ij[ind_i, ind1+1] +=(q_vec[ind1+1]-q_vec[ind2+1])\
                                   /rho0_vec[ind_i]
        ##
        DqDeta_ij[ind_i, ind2+1] = - (q_vec[ind1+1]-q_vec[ind2+1])\
                                    /rho0_vec[ind_i]

    ###############
    rho_vec = DetermineLocalCoordinates(q_vec)
    ####
    # eta_vec = np.zeros(9)
    # for ind_i in range(9):
    #     eta_vec[ind_i] = (rho_vec[ind_i]/rho0_vec[ind_i] \
    #                       - rho0_vec[ind_i])\
    #                       /(2.0)
    eta_vec = 0.5*(rho_vec/rho0_vec - rho0_vec)

    # Dq_vec = np.zeros(12)
    # for ind_i in range(12):
    #     for ind_j in range(9):
    #         for ind_k in range(9):
    #             Dq_vec[ind_i] = Dq_vec[ind_i] \
    #                             + DqDeta_ij[ind_j, ind_i] \
    #                             * M_ij[ind_j, ind_k] \
    #                             * eta_vec[ind_k]
    Dq_vec = np.dot(DqDeta_ij.T,np.dot(M_ij,eta_vec))
    ######################
    #dforce = np.zeros([6, 2])
    #for ind_q in range(6):
    #    dforce[ind_q, 0] = Dq_vec[2*ind_q]
    #    dforce[ind_q, 1] = Dq_vec[2*ind_q + 1]
    return np.array([Dq_vec[0::2],Dq_vec[1::2]]).T
########################################################
########################################################

########################################################
########################################################
def FunctionBulkEnergy(l_temp, *args):
    ###############
    (M_ij, rho0_vec)= args
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
    (M_ij, rho0_vec)= args
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
    (M_ij, rho0_vec)= args
    l0_guess = 0.9
    ###############
    res = fmin_cg(FunctionBulkEnergy, l0_guess, \
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
