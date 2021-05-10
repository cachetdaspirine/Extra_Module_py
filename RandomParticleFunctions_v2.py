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
def RandomParticle(seed_temp,pressure):
    n_t = 2
    ##########################
    np.random.seed(seed_temp)
    ## Random matrix
    m_ij = np.random.rand(12,12)
    ##########################
    for ind_i in range(12):
        m_ij[ind_i, ind_i] = m_ij[ind_i, ind_i] + n_t
    q0_vec,eps1,eps2 = RandomPositions(seed_temp)
    m_ij+=AreaMatrix(pressure)
    ##########################
    ## symmetrize
    ##########################
    ## Three-fold
    m_ij = RotateThreeFold(m_ij)
    ## i <--> j symmetry
    m_ij = MakeSymmetricMatrix(m_ij)
    ## Translational symmetry
    m_ij = ApplyTranslationalSymmetry(m_ij)
    ## Rotational symmetry
    m_ij = ApplyRotationalSymmetry(m_ij, q0_vec)
    ## Rotate m_ij and q0 to match Hugo's order of nodes
    #Exchange the index 0 into 5
    #for X
    #q0_vec = np.append(q0_vec,q0_vec[0])
    #q0_vec = np.delete(q0_vec,0)
    # for Y
    #q0_vec = np.append(q0_vec,q0_vec[0])
    #q0_vec = np.delete(q0_vec,0)
    #Exchange the column 0 into 5
    #for X
    #m_ij = np.append(m_ij,m_ij[:,0][:,np.newaxis],axis=1)
    #m_ij = np.delete(m_ij,0,axis = 1)
    #for Y
    #m_ij = np.append(m_ij,m_ij[:,0][:,np.newaxis],axis=1)
    #m_ij = np.delete(m_ij,0,axis = 1)
    #Exchange the line 0 into 5
    #for X
    #m_ij = np.append(m_ij,m_ij[0][np.newaxis,:],axis=0)
    #m_ij = np.delete(m_ij,0,axis=0)
    #for Y
    #m_ij = np.append(m_ij,m_ij[0][np.newaxis,:],axis=0)
    #m_ij = np.delete(m_ij,0,axis=0)
    ##########################
    ##########################
    return (m_ij, q0_vec,eps1,eps2)
########################################################
########################################################
########################################################
########################################################
def RandomPositions(seed_temp):
    ##########################
    np.random.seed(seed_temp)
    ## Two random epsilon
    epsilon_1 = 0.1#np.random.rand(1)[0]/10.0
    epsilon_2 = 0.#np.random.rand(1)[0]/10.0
    ## Regular Hexagon
    ell2 = 1.0 + epsilon_1
    q0_vec = np.zeros(12)
    # q0_vec[0] = 1.0 * l0
    # q0_vec[1] = 0.0
    # q0_vec[2] = ell2*math.cos(math.pi/3 + epsilon_2) * l0
    # q0_vec[3] = ell2*math.sin(math.pi/3 + epsilon_2) * l0
    # q0_vec[4] = -1.0/2.0 * l0
    # q0_vec[5] = math.sqrt(3.0/4.0) *l0
    # q0_vec[6] = ell2*math.cos(math.pi + epsilon_2) *l0
    # q0_vec[7] = ell2*math.sin(math.pi + epsilon_2)*l0
    # q0_vec[8] = -1.0/2.0*l0
    # q0_vec[9] = -math.sqrt(3.0/4.0)*l0
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
    #q0_vec = np.array([0.5774,0.,0.2887,0.5,-0.2887,0.5,-0.5774,0.,-0.2887,-0.5,0.2887,-0.5])
    return q0_vec,epsilon_1,epsilon_2
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
########################################################
def RotateThreeFold(m_ij):
    ###########
    ## Rotation matrix
    R_ij = np.zeros([12,12])
    R_ij_inv = np.zeros([12,12])
    for ind_i in range(12):
        if ind_i%2==0:
            ##
            R_ij[ind_i, ind_i] = - 1.0/2.0
            R_ij[ind_i, ind_i+1] = - math.sqrt(3.0)/2.0
            R_ij[ind_i+1, ind_i] = math.sqrt(3.0)/2.0
            R_ij[ind_i+1, ind_i+1] = - 1.0/2.0
            ##
            R_ij_inv[ind_i, ind_i] = - 1.0/2.0
            R_ij_inv[ind_i, ind_i+1] = math.sqrt(3.0)/2.0
            R_ij_inv[ind_i+1, ind_i] = - math.sqrt(3.0)/2.0
            R_ij_inv[ind_i+1, ind_i+1] = - 1.0/2.0
    ###########
    ## Permutation matrix
    P_ij = np.zeros([12,12])
    P_ij_inv = np.zeros([12,12])
    for ind_i in range(12):
        for ind_j in range(12):
            if (ind_i-4)%12==ind_j:
                P_ij[ind_i, ind_j] = 1.0
            if (ind_i+4)%12==ind_j:
                P_ij_inv[ind_i, ind_j] = 1.0
    ###########
    ## Combine both R and P
    PR_ij = np.zeros([12, 12])
    PR_ij_inv = np.zeros([12, 12])
    PR_ij = np.matmul(P_ij, R_ij)
    PR_ij_inv = np.matmul(P_ij_inv, R_ij_inv)
    ###########
    m_ij_s1 = RotateMij(PR_ij, PR_ij_inv, m_ij)
    m_ij_s2 = RotateMij(PR_ij, PR_ij_inv, m_ij_s1)
    ###########
    m_ij_s0 = np.zeros([12, 12])
    for ind_i in range(12):
        for ind_j in range(12):
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
    m_ij_s = np.zeros([12, 12])
    for ind_i in range(12):
        for ind_j in range(12):
            m_ij_s[ind_i, ind_j] = (1.0/2.0)*\
                                   (m_ij[ind_i, ind_j] \
                                    + m_ij[ind_j, ind_i])
    return m_ij_s
########################################################
########################################################
########################################################
########################################################
def ApplyTranslationalSymmetry(m_ij):
    ###########
    m_ij_s = np.zeros([12, 12])
    I_ij = np.zeros([12, 12])
    ###########
    for ind_i in range(12):
        I_ij[ind_i,ind_i] = 1.0
    ###########
    n_t = math.sqrt(1.0/6.0)
    u_x = [n_t, 0.0, n_t, 0.0, n_t, 0.0, \
           n_t, 0.0, n_t, 0.0, n_t, 0.0]
    u_y = [0.0, n_t, 0.0, n_t, 0.0, n_t, \
           0.0, n_t, 0.0, n_t, 0.0, n_t]
    ux_ij = np.outer(u_x, u_x)
    uy_ij = np.outer(u_y, u_y)
    ###########
    Trans_ij = np.zeros([12,12])
    Trans_ij_inv = np.zeros([12,12])
    for ind_i in range(12):
        for ind_j in range(12):
            Trans_ij[ind_i, ind_j] = I_ij[ind_i,ind_j] \
                                     - ux_ij[ind_i, ind_j] \
                                     - uy_ij[ind_i, ind_j]
            Trans_ij_inv[ind_j, ind_i] = I_ij[ind_i,ind_j] \
                                         - ux_ij[ind_i, ind_j] \
                                         - uy_ij[ind_i, ind_j]
    ###########
    m_ij_s0 = np.matmul(Trans_ij_inv, m_ij)
    m_ij_s = np.matmul(m_ij_s0, Trans_ij)
    ###########
    return m_ij_s
########################################################
########################################################
########################################################
########################################################
def ApplyRotationalSymmetry(m_ij, q0_vec):
    ###########
    ## Infinitesimal rotation
    ###########
    m_ij_s = np.zeros([12, 12])
    I_ij = np.zeros([12, 12])
    ###########
    for ind_i in range(12):
        I_ij[ind_i,ind_i] = 1.0
    ###########
    q0_rotated = np.zeros(12)
    norm_q0 = math.sqrt(np.dot(q0_vec, q0_vec))
    for ind_ir in range(12):
        if ind_ir%2==0:
            q0_rotated[ind_ir] = - q0_vec[ind_ir+1]\
                                 /norm_q0
        elif ind_ir%2==1:
            q0_rotated[ind_ir] = q0_vec[ind_ir-1]\
                                 /norm_q0
    ###########
    d_t_ij = np.outer(q0_rotated, q0_rotated)
    ###########
    Rot_ij = np.zeros([12,12])
    Rot_ij_inv = np.zeros([12,12])
    for ind_i in range(12):
        for ind_j in range(12):
            Rot_ij[ind_i, ind_j] = I_ij[ind_i,ind_j] \
                                   - d_t_ij[ind_i, ind_j]
            Rot_ij_inv[ind_j, ind_i] = I_ij[ind_i,ind_j] \
                                       - d_t_ij[ind_i, ind_j]
    ###########
    m_ij_s0 = np.matmul(Rot_ij_inv, m_ij)
    m_ij_s = np.matmul(m_ij_s0, Rot_ij)
    ###########
    return m_ij_s
########################################################
########################################################
########################################################
########################################################
def PrintMij(m_ij_t):
    for ind_i in range(12):
        print('%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f' \
              %(m_ij_t[ind_i, 0], m_ij_t[ind_i, 1],\
                m_ij_t[ind_i, 2], \
                m_ij_t[ind_i, 3], m_ij_t[ind_i, 4],\
                m_ij_t[ind_i, 5], \
                m_ij_t[ind_i, 6], m_ij_t[ind_i, 7],\
                m_ij_t[ind_i, 8], \
                m_ij_t[ind_i, 9], m_ij_t[ind_i, 10],\
                m_ij_t[ind_i, 11]))
########################################################
########################################################
########################################################
########################################################
def FindEigenValues(m_ij_t):
    res = np.linalg.eigh(m_ij_t)
    eign = res[0]
    #return res
    return eign
########################################################
########################################################
########################################################
########################################################
# def AreaMatrix(pressure):
    # m_ij_A = np.zeros([12,12])
    # for ind_i in range(12):
        # for ind_j in range(12):
            # temp_add = 0.0
            # if ind_i%2==0 and (ind_i+3)%12==ind_j:
                # temp_add = pressure/4.0
            # elif ind_j%2==0 and (ind_j+3)%12==ind_i:
                # temp_add = pressure/4.0
            # elif ind_i%2==1 and (ind_i+1)%12==ind_j:
                # temp_add = - pressure/4.0
            # elif ind_j%2==1 and (ind_j+1)%12==ind_i:
                # temp_add = - pressure/4.0
            # m_ij_A[ind_i, ind_j] = temp_add

def AreaMatrix(pressure):
    m_ij_A = np.zeros((12,12),dtype=float)
    for i in range(12):
        m_ij_A[i,(i+3)%12]+=(1-2*(i%2)) * pressure
    m_ij_A = (m_ij_A+np.transpose(m_ij_A))
    return m_ij_A
########################################################
########################################################
########################################################
########################################################
def AddAreaMatrices(m_ij, q0_vec, pressure):
    ##
    m_ij_A = AreaMatrix(pressure)
    ##
    m_ij_temp = np.zeros([12, 12])
    for ind_i in range(12):
        for ind_j in range(12):
            m_ij_temp[ind_i, ind_j] = m_ij[ind_i, ind_j] + \
                                      m_ij_A[ind_i, ind_j]
    ##
    A_vec_temp = m_ij_A.dot(q0_vec)
    ##
    return (m_ij_temp, A_vec_temp)
########################################################
########################################################
########################################################
########################################################
