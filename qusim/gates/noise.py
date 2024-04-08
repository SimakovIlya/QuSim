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

    II = np.sqrt(1-15*p/16)*tensor(to_Pauli_T_matrix(np.kron(I,I))).permute(l)
    IX = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(I,X))).permute(l)
    IY = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(I,Y))).permute(l)
    IZ = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(I,Z))).permute(l)
    XI = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(X,I))).permute(l)
    XX = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(X,X))).permute(l)
    XY = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(X,Y))).permute(l)
    XZ = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(X,Z))).permute(l)
    YI = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Y,I))).permute(l)
    YX = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Y,X))).permute(l)
    YY = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Y,Y))).permute(l)
    YZ = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Y,Z))).permute(l)
    ZI = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Z,I))).permute(l)
    ZX = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Z,X))).permute(l)
    ZY = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Z,Y))).permute(l)
    ZZ = np.sqrt(p/16)*tensor(to_Pauli_T_matrix(np.kron(Z,Z))).permute(l)

    return(Qobj(II+IX+IY+IZ+XI+XX+XY+XZ+YI+YX+YY+YZ+ZI+ZX+ZY+ZZ))


def R3_depolarization(p): # transfer matrix depolarizing in Pauli basis:
    # (1-p)rho + p*I/2
    I = np.array([[1,0],[0,1]], dtype = complex)
    X = np.array([[0,1],[1,0]], dtype = complex)
    Y = np.array([[0,-1j],[1j, 0]], dtype = complex)
    Z = np.array([[1,0],[0,-1]], dtype = complex)
    l = [0, 1]

    III = np.sqrt(1 - 63 * p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, I)))).permute(l)
    IIX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, X)))).permute(l)
    IIY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, Y)))).permute(l)
    IIZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(I, Z)))).permute(l)
    IXI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, I)))).permute(l)
    IXX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, X)))).permute(l)
    IXY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, Y)))).permute(l)
    IXZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(X, Z)))).permute(l)
    IYI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, I)))).permute(l)
    IYX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, X)))).permute(l)
    IYY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, Y)))).permute(l)
    IYZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Y, Z)))).permute(l)
    IZI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, I)))).permute(l)
    IZX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, X)))).permute(l)
    IZY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, Y)))).permute(l)
    IZZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(I, np.kron(Z, Z)))).permute(l)
    XII = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, I)))).permute(l)
    XIX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, X)))).permute(l)
    XIY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, Y)))).permute(l)
    XIZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(I, Z)))).permute(l)
    XXI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, I)))).permute(l)
    XXX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, X)))).permute(l)
    XXY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, Y)))).permute(l)
    XXZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(X, Z)))).permute(l)
    XYI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, I)))).permute(l)
    XYX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, X)))).permute(l)
    XYY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, Y)))).permute(l)
    XYZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Y, Z)))).permute(l)
    XZI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, I)))).permute(l)
    XZX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, X)))).permute(l)
    XZY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, Y)))).permute(l)
    XZZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(X, np.kron(Z, Z)))).permute(l)
    YII = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, I)))).permute(l)
    YIX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, X)))).permute(l)
    YIY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, Y)))).permute(l)
    YIZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(I, Z)))).permute(l)
    YXI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, I)))).permute(l)
    YXX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, X)))).permute(l)
    YXY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, Y)))).permute(l)
    YXZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(X, Z)))).permute(l)
    YYI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, I)))).permute(l)
    YYX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, X)))).permute(l)
    YYY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, Y)))).permute(l)
    YYZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Y, Z)))).permute(l)
    YZI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, I)))).permute(l)
    YZX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, X)))).permute(l)
    YZY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, Y)))).permute(l)
    YZZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Y, np.kron(Z, Z)))).permute(l)
    ZII = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, I)))).permute(l)
    ZIX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, X)))).permute(l)
    ZIY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, Y)))).permute(l)
    ZIZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(I, Z)))).permute(l)
    ZXI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, I)))).permute(l)
    ZXX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, X)))).permute(l)
    ZXY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, Y)))).permute(l)
    ZXZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(X, Z)))).permute(l)
    ZYI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, I)))).permute(l)
    ZYX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, X)))).permute(l)
    ZYY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, Y)))).permute(l)
    ZYZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Y, Z)))).permute(l)
    ZZI = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, I)))).permute(l)
    ZZX = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, X)))).permute(l)
    ZZY = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, Y)))).permute(l)
    ZZZ = np.sqrt(p / 64) * tensor(to_Pauli_T_matrix(np.kron(Z, np.kron(Z, Z)))).permute(l)

    return(Qobj(III+IIX+IIY+IIZ+IXI+IXX+IXY+IXZ+IYI+IYX+IYY+IYZ+IZI+IZX+IZY+IZZ+
                XII+XIX+XIY+XIZ+XXI+XXX+XXY+XXZ+XYI+XYX+XYY+XYZ+XZI+XZX+XZY+XZZ+
                YII+YIX+YIY+YIZ+YXI+YXX+YXY+YXZ+YYI+YYX+YYY+YYZ+YZI+YZX+YZY+YZZ+
                ZII+ZIX+ZIY+ZIZ+ZXI+ZXX+ZXY+ZXZ+ZYI+ZYX+ZYY+ZYZ+ZZI+ZZX+ZZY+ZZZ))



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








