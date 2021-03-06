from ..Pauli_basis import *

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




def iSWAP_gate(theta = np.pi, eta = 0, theta_error = 0, eta_error = 0, phi_error = 0):
    '''
    Returns Pauli trasfer matrix of iSWAP gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, np.cos(theta/2)*np.exp(-theta_error**2/8)*np.exp(-phi_error**2/2), 1j*np.exp(1j*eta)*np.exp(-eta_error**2/2)*np.exp(-theta_error**2/8)*np.sin(theta/2), 0],
                                    [0, 1j*np.exp(-1j*eta)*np.exp(-eta_error**2/2)*np.sin(theta/2)*np.exp(-theta_error**2/8), np.cos(theta/2)*np.exp(-theta_error**2/8)*np.exp(-phi_error**2/2), 0],
                                    [0, 0, 0, np.exp(-phi_error**2/8)]]))
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