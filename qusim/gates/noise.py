from ..Pauli_basis import *


def R_dep(p_plane, p_axis): # transfer matrix depolarizing in Pauli basis
    return(Qobj(np.array([[1, 0, 0, 0],
                          [0, 1 - p_plane, 0, 0],
                          [0, 0, 1 - p_axis, 0],
                          [0, 0, 0, 1 - p_plane]])))


def R_dephasing_2qubits_mse(phi): # transfer dephasing matrix for 2 qubits in Pauli basis with mean sq. phase error
    a = np.exp(-phi**2/2)
    A = Qobj(np.identity(4))
    A = tensor(A, A)
    for i in range(4, 12):
        A.data[i, i] = a
    return(A)


def R_dephasing_11_mse(phi): # transfer dephasing matrix for 2 qubits in Pauli basis with mean square phase error
    a = np.exp(-phi**2/2)
    A = Qobj(np.identity(4))
    A = tensor(A, A)
    A.data[3, 3] = a
    A.data[7, 7] = a
    A.data[11, 11] = a
    A.data[12, 12] = a
    A.data[13, 13] = a
    A.data[14, 14] = a
    return(A)


def add_noise_Ry(Ry, dtheta = 0):
    '''
    PTM Ry(theta) (Qobj):
    1,      0,      0,     0
    0, cos(theta),  0, sin(theta)
    0,      0,      1,     0
    0, -sin(theta), 0, cos(theta)
    '''
    a = np.exp(-dtheta**2/2)
    Ry.data[1, 1] *= a
    Ry.data[1, 3] *= a
    Ry.data[3, 1] *= a
    Ry.data[3, 3] *= a
    return(Ry)



def add_noise_Rx(Rx, dtheta = 0):
    '''
    PTM Rx(theta) (Qobj):
    1,      0,      0,     0
    0,      0,      0,     0
    0,      0, cos(theta), sin(theta)
    0,      0, -sin(theta),cos(theta)
    '''
    a = np.exp(-dtheta**2/2)
    Rx.data[2, 2] *= a
    Rx.data[2, 3] *= a
    Rx.data[3, 2] *= a
    Rx.data[3, 3] *= a
    return(Rx)





def add_noise_iSWAP(gate, dtheta = 0, deta = 0):

    a = np.exp(-deta**2/2) * np.exp(-dtheta**2/8)
    b = 1/2*(np.exp(-dtheta**2/2) + 1)
    c = 1/2 + 1/4*(np.exp(-dtheta**2/2)) + 1/4*(np.exp(-deta**2/2))

    gate.data[1, 11]  *= a
    gate.data[2, 7]   *= a
    gate.data[3, 12]  *= b
    gate.data[4, 14]  *= a
    gate.data[5, 5]   *= c
    gate.data[6, 9]   *= c
    gate.data[7, 2]   *= a
    gate.data[8, 13]  *= a
    gate.data[9, 6]   *= c
    gate.data[10, 10] *= c
    gate.data[11, 1]  *= a
    gate.data[12, 4]  *= b
    gate.data[13, 8]  *= a
    gate.data[14, 4]  *= a

    return(gate)



def add_noise_ZZ(gate, dphi = 0):

    a = 1/2*np.exp(-dphi**2/2) + 1/2*np.exp(-dphi**2/8)
    b = 1/2*(np.exp(-dphi**2/8) + 1)

    gate.data[1, 1]   *= a
    gate.data[2, 2]   *= a
    gate.data[4, 4]   *= a
    gate.data[5, 5]   *= b
    gate.data[6, 6]   *= b
    gate.data[7, 7]   *= a
    gate.data[8, 8]   *= a
    gate.data[9, 9]   *= b
    gate.data[10, 10] *= b
    gate.data[11, 11] *= a
    gate.data[13, 13] *= a
    gate.data[14, 14] *= a

    return(gate)








