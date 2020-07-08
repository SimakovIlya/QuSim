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


def add_noise_Ry(Ry, dtheta):
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