from .basis import *




def to_Pauli_T_matrix(O):
    '''
    Converts one or two qubit gate matrix to Pauli transfer matrix as Qobj
    
    Args:
        ndarray: gate in ordinary basis (2x2 or 4x4)
        
    Returns:
        ndarray: Pauli transfer matrix (4x4 or 16x16)
    '''
    if (O.shape[0] == 2):
        T = Qobj(np.eye(4))
        for i in range(0, 4):
            for j in range(0, 4):
                T.data[i, j] = 1/2 * np.trace(sigma[i]@O@sigma[j]@np.conj(O.T))
    if (O.shape[0] == 4):
        T = tensor(Qobj(np.eye(4)), Qobj(np.eye(4)))
        sigma2d = np.empty((16, 4, 4), np.complex)
        for i in range(0, 4):
            for j in range(0, 4):
                sigma2d[4*i + j,:,:] = np.kron(sigma[i], sigma[j])
        for i in range(0, 16):
            for j in range(0, 16):
                T.data[i, j] = 1/4 * np.trace(sigma2d[i,:,:]@O@sigma2d[j,:,:]@np.transpose(O.conj()))       
    return(T)