from import_Pauli_basis import *


def CZ_gate():
    '''
    Returns trasfer matrix of Ry gate
    '''
    l = [0, 1]
    T = to_Pauli_T_matrix(np.array([[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, -1]]))
    return tensor(T).permute(l)