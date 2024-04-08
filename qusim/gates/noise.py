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



def R2_depolarization(p): # transfer matrix depolarizing in Pauli basis:
    # (1-p)rho + p*I/2
    I = np.array([[1,0],[0,1]], dtype = complex)
    X = np.array([[0,1],[1,0]], dtype = complex)
    Y = np.array([[0,-1j],[1j, 0]], dtype = complex)
    Z = np.array([[1,0],[0,-1]], dtype = complex)
    l = [0, 1]

    # II = (1 - 15 * p / 16) * tensor(to_Pauli_T_matrix(np.kron(I, I))).permute(l)
    # IX = p / 16 * tensor(to_Pauli_T_matrix(np.kron(I, X))).permute(l)
    # IY = p / 16 * tensor(to_Pauli_T_matrix(np.kron(I, Y))).permute(l)
    # IZ = p / 16 * tensor(to_Pauli_T_matrix(np.kron(I, Z))).permute(l)
    # XI = p / 16 * tensor(to_Pauli_T_matrix(np.kron(X, I))).permute(l)
    # XX = p / 16 * tensor(to_Pauli_T_matrix(np.kron(X, X))).permute(l)
    # XY = p / 16 * tensor(to_Pauli_T_matrix(np.kron(X, Y))).permute(l)
    # XZ = p / 16 * tensor(to_Pauli_T_matrix(np.kron(X, Z))).permute(l)
    # YI = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Y, I))).permute(l)
    # YX = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Y, X))).permute(l)
    # YY = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Y, Y))).permute(l)
    # YZ = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Y, Z))).permute(l)
    # ZI = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Z, I))).permute(l)
    # ZX = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Z, X))).permute(l)
    # ZY = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Z, Y))).permute(l)
    # ZZ = p / 16 * tensor(to_Pauli_T_matrix(np.kron(Z, Z))).permute(l)
    #
    # return(Qobj(II+IX+IY+IZ+XI+XX+XY+XZ+YI+YX+YY+YZ+ZI+ZX+ZY+ZZ))

    R = tensor(to_Pauli_T_matrix(np.kron(I, I))).permute(l)
    for i in range(1, R.shape[0]):
        R.data[i,i] = 1-p
    return R


def R3_depolarization(p): # transfer matrix depolarizing in Pauli basis:
    # (1-p)rho + p*I/2
    I = np.array([[1,0],[0,1]], dtype = complex)
    X = np.array([[0,1],[1,0]], dtype = complex)
    Y = np.array([[0,-1j],[1j, 0]], dtype = complex)
    Z = np.array([[1,0],[0,-1]], dtype = complex)
    l = [0, 1, 2]

    # III = (1 - 63 * p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, I)))).permute(l)
    # IIX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, X)))).permute(l)
    # IIY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, Y)))).permute(l)
    # IIZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, Z)))).permute(l)
    # IXI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, I)))).permute(l)
    # IXX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, X)))).permute(l)
    # IXY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, Y)))).permute(l)
    # IXZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, Z)))).permute(l)
    # IYI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, I)))).permute(l)
    # IYX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, X)))).permute(l)
    # IYY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, Y)))).permute(l)
    # IYZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, Z)))).permute(l)
    # IZI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, I)))).permute(l)
    # IZX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, X)))).permute(l)
    # IZY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, Y)))).permute(l)
    # IZZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, Z)))).permute(l)
    # XII = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, I)))).permute(l)
    # XIX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, X)))).permute(l)
    # XIY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, Y)))).permute(l)
    # XIZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, Z)))).permute(l)
    # XXI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, I)))).permute(l)
    # XXX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, X)))).permute(l)
    # XXY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, Y)))).permute(l)
    # XXZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, Z)))).permute(l)
    # XYI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, I)))).permute(l)
    # XYX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, X)))).permute(l)
    # XYY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, Y)))).permute(l)
    # XYZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, Z)))).permute(l)
    # XZI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, I)))).permute(l)
    # XZX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, X)))).permute(l)
    # XZY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, Y)))).permute(l)
    # XZZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, Z)))).permute(l)
    # YII = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, I)))).permute(l)
    # YIX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, X)))).permute(l)
    # YIY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, Y)))).permute(l)
    # YIZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, Z)))).permute(l)
    # YXI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, I)))).permute(l)
    # YXX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, X)))).permute(l)
    # YXY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, Y)))).permute(l)
    # YXZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, Z)))).permute(l)
    # YYI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, I)))).permute(l)
    # YYX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, X)))).permute(l)
    # YYY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, Y)))).permute(l)
    # YYZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, Z)))).permute(l)
    # YZI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, I)))).permute(l)
    # YZX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, X)))).permute(l)
    # YZY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, Y)))).permute(l)
    # YZZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, Z)))).permute(l)
    # ZII = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, I)))).permute(l)
    # ZIX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, X)))).permute(l)
    # ZIY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, Y)))).permute(l)
    # ZIZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, Z)))).permute(l)
    # ZXI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, I)))).permute(l)
    # ZXX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, X)))).permute(l)
    # ZXY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, Y)))).permute(l)
    # ZXZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, Z)))).permute(l)
    # ZYI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, I)))).permute(l)
    # ZYX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, X)))).permute(l)
    # ZYY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, Y)))).permute(l)
    # ZYZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, Z)))).permute(l)
    # ZZI = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, I)))).permute(l)
    # ZZX = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, X)))).permute(l)
    # ZZY = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, Y)))).permute(l)
    # ZZZ = p / 64 * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, Z)))).permute(l)
    #
    # return(Qobj(III+IIX+IIY+IIZ+IXI+IXX+IXY+IXZ+IYI+IYX+IYY+IYZ+IZI+IZX+IZY+IZZ+
    #             XII+XIX+XIY+XIZ+XXI+XXX+XXY+XXZ+XYI+XYX+XYY+XYZ+XZI+XZX+XZY+XZZ+
    #             YII+YIX+YIY+YIZ+YXI+YXX+YXY+YXZ+YYI+YYX+YYY+YYZ+YZI+YZX+YZY+YZZ+
    #             ZII+ZIX+ZIY+ZIZ+ZXI+ZXX+ZXY+ZXZ+ZYI+ZYX+ZYY+ZYZ+ZZI+ZZX+ZZY+ZZZ))
    R = tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, I)))).permute(l)
    for i in range(1, R.shape[0]):
        R.data[i,i] = 1-p
    return R



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








