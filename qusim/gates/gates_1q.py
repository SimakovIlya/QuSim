from ..Pauli_basis import *
    

def SingleQubit_gate(U):
    '''
    Returns trasfer matrix of U gate
    '''
    return to_Pauli_T_matrix(np.asarray(U))




def NOT_gate():
    '''
    Returns Pauli trasfer matrix of NOT gate
    '''
    return to_Pauli_T_matrix(np.array([[0, 1],
                                       [1, 0]]))



def Ry_gate(theta):
    '''
    Returns Pauli trasfer matrix of Ry(theta) gate
    '''
    return to_Pauli_T_matrix(np.array([[np.cos(theta/2), -np.sin(theta/2)],
                                       [np.sin(theta/2), np.cos(theta/2)]]))



def Rx_gate(theta):
    '''
    Returns Pauli trasfer matrix of Rx(theta) gate
    '''
    return to_Pauli_T_matrix(np.array([[np.cos(theta/2), -1j*np.sin(theta/2)],
                                       [-1j*np.sin(theta/2), np.cos(theta/2)]]))




def Rz_gate(theta):
    '''
    Returns Pauli trasfer matrix of Rx(theta) gate
    '''
    return to_Pauli_T_matrix(np.array([[1, 0],
                                       [0, np.exp(1j*theta)]]))




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



def X_gate():
    '''
    Returns Pauli trasfer matrix of X gate
    '''
    return to_Pauli_T_matrix(np.array([[0, 1],
                                       [1, 0]]))




def Y_gate():
    '''
    Returns Pauli trasfer matrix of Y gate
    '''
    return to_Pauli_T_matrix(np.array([[0, -1j],
                                       [1j, 0]]))