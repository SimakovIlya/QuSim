from ..Pauli_basis import *
    

def NOT_gate():
    '''
    Returns Pauli trasfer matrix of NOT gate
    '''
    return to_Pauli_T_matrix(np.array([[0, 1],
                                       [1, 0]]))



def Ry_gate(theta, theta_error = 0):
    '''
    Returns Pauli trasfer matrix of Ry(theta) gate
    '''
    return to_Pauli_T_matrix(np.array([[np.cos(theta/2)*np.exp(-1j*theta_error**2/8), -np.sin(theta/2)],
                                       [np.sin(theta/2), np.cos(theta/2)*np.exp(-1j*theta_error**2/8)]]))
    # Qobj(np.array([[1, 0, 0, 0],
    #                       [0, np.cos(theta)*np.exp(-1j*theta_error**2/8), 0, np.sin(theta)],
    #                       [0, 0, 1, 0],
    #                       [0, -np.sin(theta), 0, np.cos(theta)*np.exp(-1j*theta_error**2/8)]]))




def Hadamard_gate():
    '''
    Returns Pauli trasfer matrix of Hadamard gate
    '''
    return to_Pauli_T_matrix(1/np.sqrt(2)*np.array([[1, 1],
                                                    [1, -1]]))




def S_gate(sign = 1):
    '''
    Returns Pauli trasfer matrix of S phase gate
    '''
    return to_Pauli_T_matrix(np.array([[1, 0],
                                       [0, sign * 1j]]))



def Z_gate():
    '''
    Returns Pauli trasfer matrix of Pauli-Z gate
    '''
    return to_Pauli_T_matrix(np.array([[1, 0],
                                       [0, -1]]))