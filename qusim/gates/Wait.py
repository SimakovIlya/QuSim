from ..Pauli_basis import *

def Wait(t, T1, Tf):
    '''
    Returns trasfer matrix of idling
    
    Args:
        t (float): time of the idling 
        T1 (float): T1
        Tf (float): T_phi
        
    Returns:
        Qobj: trasfer matrix of idling (4 x 4)
    '''
    R = Qobj(RT1(t, T1)@RTf(t, Tf))
    return(R)




def RT1(t, T1): # transfer matrix amplitude damping in Pauli basis
    #print(t, T1)
    p = 1 - np.exp(-t/T1)
    return(np.array([[1, 0, 0, 0],
                     [0, np.sqrt(1 - p), 0, 0],
                     [0, 0, np.sqrt(1 - p), 0],
                     [p, 0, 0, 1 - p]]))


def RTf(t, Tf): # transfer matrix phase damping in Pauli basis
    p = 1 - np.exp(-2*t/Tf)
    return(np.array([[1, 0, 0, 0],
                     [0, np.sqrt(1 - p), 0, 0],
                     [0, 0, np.sqrt(1 - p), 0],
                     [0, 0, 0, 1]]))