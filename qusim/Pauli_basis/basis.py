import numpy as np
from qutip import *
from scipy.sparse import *




sigma = np.array([[[1, 0], [0, 1]],
                  [[0, 1], [1, 0]],
                  [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype = complex)


T_matrix = np.array([[1, 0, 0, 1],
                     [0, 1, -1j, 0],
                     [0, 1, 1j, 0],
                     [1, 0, 0, -1]]) # transformation matrix from odinary basis to Pauli basis
                                     # vec(in ordinary) = trans_matrix * vec (in Pauli)
T_matrix_inv = np.linalg.inv(T_matrix)




def vector_to_Pauli(rho):
    '''
    Returns rho-vector in Pauli basis 
    
    Args:
        Qobj: 1-2dim rho-vector in ordinary basis
        
    Returns:
        Qobj: 1-2dim rho-vector in Pauli basis
        
    '''
    TM1_inv = Qobj(T_matrix_inv)
    T2_inv = tensor([TM1_inv] * 2)
    if rho.shape[0] == 4:
        return TM1_inv * rho
    else:
        return T2_inv * rho




def vector_from_Pauli(rhoP):
    '''
    Returns rho-vector in ordinary basis 
    
    Args:
        Qobj: 1dim rho-vector in Pauli basis
        
    Returns:
        Qobj: 1-2dim rho-vector in ordinary basis
    '''
    TM1 = Qobj(T_matrix)
    T2 = tensor([TM1] * 2)
    if rhoP.shape[0] == 4:
        return TM1 * rhoP
    else:
        return T2 * rhoP