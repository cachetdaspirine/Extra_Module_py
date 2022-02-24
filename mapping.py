import numpy as np
def UseMapping():
    mapping = np.zeros([2,9], dtype=int)
    P_map_ij = np.zeros([9,9])
    P_map_ij_inv = np.zeros([9,9])
    m_temp = 2
    if m_temp == 1:
        ## 1
        mapping = np.array([[0, 2, 4, 6, 8, 10, 0, 2, 4],\
                            [2, 4, 6, 8, 10, 0, 6, 8, 10]])
    elif m_temp == 2:
        ## 2
        mapping = np.array([[0, 2, 4, 6, 8,10, 0, 4, 8],\
                            [4, 6, 8,10, 0, 2, 2, 6,10]])
    elif m_temp == 3:
        ## 3
        mapping = np.array([[0, 2, 4, 6, 8,10, 0, 4, 2],\
                            [4, 6, 8,10, 0, 2, 6,10, 8]])
    return mapping
