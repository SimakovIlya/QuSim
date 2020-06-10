from ..Pauli_basis import *



def average_Pauli_gate(PTM_matrix, args):
    '''
    Calculates average Pauli gate for gate.
    
    Args:
        PTM_matrix (Qobj): Pauli Transfer Matrix of the gate
        args (list): list of ndarrays (or lists) of params for PTM, dim = [averaging number, number of params]
    Returns:
        (Qobj) Pauli Transfer Matrix of the average gate
    '''
    args = np.asarray(args)
    n = len(args[0])
    args_tmp = args[:,0]
    result = PTM_matrix(*args_tmp)
    for i in range(1, n):
        args_tmp = args[:,i]
        result +=  PTM_matrix(*args_tmp)
    return result/n