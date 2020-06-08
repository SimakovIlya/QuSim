from import_Pauli_basis import *


def CZ_gate(phi = np.pi):
    '''
    Returns trasfer matrix of Ry gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, np.exp(-1j*phi)]]))
    return tensor(T).permute(l)




def iSWAP_gate():
    '''
    Returns Pauli trasfer matrix of iSWAP gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, 0, 1j, 0],
                                    [0, 1j, 0, 0],
                                    [0, 0, 0, 1]]))
    return tensor(T).permute(l)




def CNOT_gate():
    '''
    Returns Pauli trasfer matrix of iSWAP gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 0, 1],
                                    [0, 0, 1, 0]]))
    return tensor(T).permute(l)