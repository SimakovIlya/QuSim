from ..Pauli_basis import *


def R_dep(p_plane, p_axis): # transfer matrix depolarizing in Pauli basis
    return(Qobj(np.array([[1, 0, 0, 0],
                          [0, 1 - p_plane, 0, 0],
                          [0, 0, 1 - p_axis, 0],
                          [0, 0, 0, 1 - p_plane]])))



def R_depolarization(p): # transfer matrix depolarizing in Pauli basis: (1-p)rho + p*I/2
    I = to_Pauli_T_matrix(np.sqrt(1-3*p/4)*np.array([[1, 0],
                                                    [0, 1]]))
    X = to_Pauli_T_matrix(np.sqrt(p/4)*np.array([[0, 1],
                                                 [1, 0]]))
    Y = to_Pauli_T_matrix(np.sqrt(p/4)*np.array([[0, -1j],
                                                 [1j, 0]]))
    Z = to_Pauli_T_matrix(np.sqrt(p/4)*np.array([[1, 0],
                                                 [0, -1]]))
    return(Qobj(I+X+Y+Z))


def R_depolarization_Pauli_equal(p): # transfer matrix depolarizing in Pauli basis: (1-p)I + p/3*X + p/3*Y + p/3*Z
    I = to_Pauli_T_matrix(np.sqrt(1-p)*np.array([[1, 0],
                                                    [0, 1]]))
    X = to_Pauli_T_matrix(np.sqrt(p/3)*np.array([[0, 1],
                                                 [1, 0]]))
    Y = to_Pauli_T_matrix(np.sqrt(p/3)*np.array([[0, -1j],
                                                 [1j, 0]]))
    Z = to_Pauli_T_matrix(np.sqrt(p/3)*np.array([[1, 0],
                                                 [0, -1]]))
    return(Qobj(I+X+Y+Z))

def R_depolarization_Pauli(px, py, pz): # transfer matrix depolarizing in Pauli basis: 
    # (1-px-py-pz)I + px*X + py*Y + pz*Z
    I = to_Pauli_T_matrix(np.sqrt(1-px-py-pz)*np.array([[1, 0],
                                                    [0, 1]]))
    X = to_Pauli_T_matrix(np.sqrt(px/15)*np.array([[0, 1],
                                                 [1, 0]]))
    Y = to_Pauli_T_matrix(np.sqrt(py/15)*np.array([[0, -1j],
                                                 [1j, 0]]))
    Z = to_Pauli_T_matrix(np.sqrt(pz/15)*np.array([[1, 0],
                                                 [0, -1]]))
    return(Qobj(I+X+Y+Z))

def R2_depolarization_Pauli(p): # transfer matrix depolarizing in Pauli basis: 
    # (1-p)I + p/15*IX + p/15*IY +...+ p/15*ZZ
    I = np.array([[1,0],[0,1]], dtype = complex)
    X = np.array([[0,1],[1,0]], dtype = complex)
    Y = np.array([[0,-1j],[1j, 0]], dtype = complex)
    Z = np.array([[1,0],[0,-1]], dtype = complex)
    l = [0, 1]

    II = np.sqrt(1-p)*tensor(to_Pauli_T_matrix(np.kron(I,I))).permute(l)
    IX = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(I,X))).permute(l)
    IY = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(I,Y))).permute(l)
    IZ = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(I,Z))).permute(l)
    XI = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(X,I))).permute(l)
    XX = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(X,X))).permute(l)
    XY = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(X,Y))).permute(l)
    XZ = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(X,Z))).permute(l)
    YI = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Y,I))).permute(l)
    YX = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Y,X))).permute(l)
    YY = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Y,Y))).permute(l)
    YZ = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Y,Z))).permute(l)
    ZI = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Z,I))).permute(l)
    ZX = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Z,X))).permute(l)
    ZY = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Z,Y))).permute(l)
    ZZ = np.sqrt(p/15)*tensor(to_Pauli_T_matrix(np.kron(Z,Z))).permute(l)

    return(Qobj(II+IX+IY+IZ+XI+XX+XY+XZ+YI+YX+YY+YZ+ZI+ZX+ZY+ZZ))



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
    1,          0,           sin(theta),  0
    0,          cos(theta),  0,           0
    sin(theta), 0,           1,           0
    0,          0,           0,           cos(theta)
    '''
    a = np.exp(-dtheta**2/2)
    Rx.data[0, 3] *= a
    Rx.data[1, 1] *= a
    Rx.data[2, 0] *= a
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








