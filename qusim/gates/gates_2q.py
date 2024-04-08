from ..Pauli_basis import *

def TwoQubit_gate(U):
    '''
    Returns trasfer matrix of U gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.asarray(U))
    return tensor(T).permute(l)


def MultiQubit_gate(U):
    '''
    Returns trasfer matrix of U gate
    '''
    return TwoQubit_gate(U)

    


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




def iSWAP_gate(theta = np.pi, eta = 0):
    '''
    Returns Pauli trasfer matrix of iSWAP gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, np.cos(theta/2), 1j*np.exp(1j*eta)*np.sin(theta/2), 0],
                                    [0, 1j*np.exp(-1j*eta)*np.sin(theta/2), np.cos(theta/2), 0],
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




def fSim_gate(theta, phi):
    '''
    Returns Pauli trasfer matrix of fSim gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, np.cos(theta), -1j*np.sin(theta), 0],
                                    [0, -1j*np.sin(theta), np.cos(theta), 0],
                                    [0, 0, 0, np.exp(-1j*phi)]]))
    return tensor(T).permute(l)