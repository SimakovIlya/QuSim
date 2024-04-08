from .basis import *


def to_Pauli_T_matrix(O):
    '''
    Converts one or two qubit gate matrix to Pauli transfer matrix as ndarray

    Args:
        ndarray: gate in ordinary basis (2^n x 2^n) n=1,2,3,4

    Returns:
        ndarray: Pauli transfer matrix (4^n x 4^n)
    '''
    if (O.shape[0] == 2):
        T = np.identity(4, dtype=complex)
        for i in range(0, 4):
            for j in range(0, 4):
                T[i, j] = 1 / 2 * np.trace(sigma[i] @ O @ sigma[j] @ np.conj(O.T))
    if (O.shape[0] == 4):
        T = np.identity(16, dtype=complex)
        sigma2d = np.empty((16, 4, 4), complex)
        for i in range(0, 4):
            for j in range(0, 4):
                sigma2d[4 * i + j, :, :] = np.kron(sigma[i], sigma[j])
        for i in range(0, 16):
            for j in range(0, 16):
                T[i, j] = 1 / 4 * np.trace(sigma2d[i, :, :] @ O @ sigma2d[j, :, :] @ np.transpose(O.conj()))
    if (O.shape[0] == 8):
        T = np.identity(64, dtype=complex)
        for i in range(0, 64):
            for j in range(0, 64):
                T[i, j] = 1 / 8 * np.trace(sigma3d[i, :, :] @ O @ sigma3d[j, :, :] @ np.transpose(O.conj()))
    if (O.shape[0] == 16):
        T = np.identity(256, dtype=complex)
        for i in range(0, 256):
            for j in range(0, 256):
                T[i, j] = 1 / 16 * np.trace(sigma4d[i, :, :] @ O @ sigma4d[j, :, :] @ np.transpose(O.conj()))
    return np.real(T)


sigma2d = np.empty((16, 4, 4), complex)
for i in range(0, 4):
    for j in range(0, 4):
        sigma2d[4 * i + j, :, :] = np.kron(sigma[i], sigma[j])

sigma3d = np.empty((4 ** 3, 8, 8), complex)
for i in range(4):
    for j in range(4):
        for k in range(4):
            sigma3d[16 * i + 4 * j + k] = np.kron(np.kron(sigma[i], sigma[j]), sigma[k])

sigma4d = np.empty((4 ** 4, 16, 16), complex)
for i in range(4):
    for j in range(4):
        for k in range(4):
            for m in range(4):
                sigma4d[64 * i + 16 * j + 4 * k + m] = np.kron(np.kron(np.kron(sigma[i], sigma[j]), sigma[k]), sigma[m])