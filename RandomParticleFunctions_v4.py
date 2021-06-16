########################################################
## Functions that are used for the comparison between  #
## particle model and continuum model                  #
########################################################
########################################################
import math
from scipy.optimize import minimize
import numpy as np
########################################################
########################################################
########################################################
########################################################
########################################################
def RandomParticle(seed_temp, n_t):
    ##########################
    np.random.seed(seed_temp)
    ## Random matrix
    m_ij = np.random.rand(9,9)
    ##########################
    for ind_i in range(9):
        m_ij[ind_i, ind_i] = m_ij[ind_i, ind_i] + n_t
    rho0_vec = RandomPositions(seed_temp)
    ##########################
    ## symmetrize
    ##########################
    ## Three-fold
    m_ij = RotateThreeFold(m_ij)
    ## i <--> j symmetry
    m_ij = MakeSymmetricMatrix(m_ij)
    ##########################
    ##########################
    return (m_ij, rho0_vec)
########################################################
########################################################
########################################################
########################################################
def RandomPositions(seed_temp):
    ##########################
    np.random.seed(seed_temp)
    ## Two random epsilon
    epsilon_1 = 0.01#np.random.rand(1)/10.0
    epsilon_2 = 0.#np.random.rand(1)/10.0
    ## Regular Hexagon
    ell2 = 1.0 + epsilon_1
    q0_vec = np.zeros(12)

    q0_vec[0] = 0.5774 * (1+epsilon_1) * math.cos(math.pi/3+epsilon_2)
    q0_vec[1] = 0.5774 * (1+epsilon_1) * math.sin(math.pi/3+epsilon_2)

    q0_vec[2] = - 0.5774 * (1-epsilon_1) * math.cos(math.pi/3-epsilon_2)
    q0_vec[3] = 0.5774 * (1-epsilon_1) * math.sin(math.pi/3-epsilon_2)

    q0_vec[4] = 0.5774*(1+epsilon_1)*math.cos(math.pi+epsilon_2)
    q0_vec[5] = 0.5774*(1+epsilon_1)*math.sin(math.pi+epsilon_2)

    q0_vec[6] = -0.5774 * (1-epsilon_1) * math.cos(math.pi/3-epsilon_2)
    q0_vec[7] = -0.5774 * (1-epsilon_1) * math.sin(math.pi/3-epsilon_2)

    q0_vec[8] = 0.5774 * (1+epsilon_1) * math.cos(math.pi/3+epsilon_2)
    q0_vec[9] = -0.5774 * (1+epsilon_1) * math.sin(math.pi/3+epsilon_2)

    q0_vec[10] = 0.5774 * math.cos(2*math.pi-epsilon_2) * (1-epsilon_1)
    q0_vec[11] = 0.5774 * math.sin(2*math.pi-epsilon_2) * (1-epsilon_1)

    ##
    #mapping = np.zeros([2,9], dtype=int)
    #mapping = np.array([[0, 2, 4, 6, 8, 10, 0, 2, 4],\
    #                 [2, 4, 6, 8, 10, 0, 6, 8, 10]])
    mapping = np.array([[0,4,8,0,4,8,2,6,10,]
                        ,[2,6,10,4,8,0,6,10,2]])
    ##
    rho0_vec = np.zeros(9)
    for ind_i in range(9):
        ind1 = mapping[0, ind_i]
        ind2 = mapping[1, ind_i]
        rho0_vec[ind_i] = math.sqrt((q0_vec[ind1]-q0_vec[ind2])\
                                    *(q0_vec[ind1]-q0_vec[ind2])\
                                    +(q0_vec[ind1+1]-q0_vec[ind2+1])\
                                    *(q0_vec[ind1+1]-q0_vec[ind2+1]))
    ##
    return rho0_vec
########################################################
########################################################
########################################################
########################################################
def RotateThreeFold(m_ij):
    ###########
    ###########
    ## Permutation matrix
    P_ij = np.zeros([9, 9])
    P_ij[0, 4] = 1.0
    P_ij[1, 5] = 1.0
    P_ij[2, 0] = 1.0
    P_ij[3, 1] = 1.0
    P_ij[4, 2] = 1.0
    P_ij[5, 3] = 1.0
    P_ij[6, 8] = 1.0
    P_ij[7, 6] = 1.0
    P_ij[8, 7] = 1.0
    P_ij_inv = np.zeros([9,9])
    P_ij_inv[0, 2] = 1.0
    P_ij_inv[1, 3] = 1.0
    P_ij_inv[2, 4] = 1.0
    P_ij_inv[3, 5] = 1.0
    P_ij_inv[4, 0] = 1.0
    P_ij_inv[5, 1] = 1.0
    P_ij_inv[6, 7] = 1.0
    P_ij_inv[7, 8] = 1.0
    P_ij_inv[8, 6] = 1.0
    ###########
    m_ij_s1 = RotateMij(P_ij, P_ij_inv, m_ij)
    m_ij_s2 = RotateMij(P_ij, P_ij_inv, m_ij_s1)
    ###########
    m_ij_s0 = np.zeros([9, 9])
    for ind_i in range(9):
        for ind_j in range(9):
            m_ij_s0[ind_i, ind_j] = (1.0/3.0)*\
                                    (m_ij[ind_i, ind_j] \
                                     + m_ij_s1[ind_i, ind_j] \
                                     + m_ij_s2[ind_i, ind_j])
    ##########
    return m_ij_s0
########################################################
########################################################
########################################################
########################################################
def MakeSymmetricMatrix(m_ij):
    ###########
    ## The i <--> j symmetry
    m_ij_s = np.zeros([9, 9])
    for ind_i in range(9):
        for ind_j in range(9):
            m_ij_s[ind_i, ind_j] = (1.0/2.0)*\
                                   (m_ij[ind_i, ind_j] \
                                    + m_ij[ind_j, ind_i])
    return m_ij_s
########################################################
########################################################
########################################################
########################################################
def RotateMij(R_ij, R_ij_inv, m_ij):
    m_ij_rot_temp = np.matmul(R_ij, m_ij)
    m_ij_rot = np.matmul(m_ij_rot_temp, R_ij_inv)
    return m_ij_rot
########################################################
########################################################
########################################################
def PrintMij(m_ij_t):
    for ind_i in range(9):
        print('%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f' \
              %(m_ij_t[ind_i, 0], m_ij_t[ind_i, 1],\
                m_ij_t[ind_i, 2], \
                m_ij_t[ind_i, 3], m_ij_t[ind_i, 4],\
                m_ij_t[ind_i, 5], \
                m_ij_t[ind_i, 6], m_ij_t[ind_i, 7],\
                m_ij_t[ind_i, 8]))
########################################################
########################################################
########################################################
########################################################
def FindEigenValues(m_ij_t):
    res = np.linalg.eigh(m_ij_t)
    eign = res[0]
    return eign
########################################################
########################################################
########################################################
########################################################
def AddTwoMatrices(m1_ij, m2_ij, pressure):
    m_ij = np.zeros([9, 9])
    for ind_i in range(9):
        for ind_j in range(9):
            m_ij[ind_i, ind_j] = m1_ij[ind_i, ind_j] \
                                 - pressure\
                                 *m2_ij[ind_i, ind_j]
    return m_ij
########################################################
########################################################
########################################################
########################################################
