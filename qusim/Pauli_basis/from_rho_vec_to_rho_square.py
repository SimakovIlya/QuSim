from .basis import *

def from_rho_vec_to_rho_sq(rho_vec, mask_to_sq_rho):
    '''
    Returns rho-matrix in ordinary basis 
    
    Args:
        Qobj: rho-vector in ordinary basis
        
    Returns:
        Qobj: rho-matrix in ordinary basis 
    '''
    n = int(np.log2(rho_vec.shape[0])/2)
    rho = Qobj(rho_vec.full().ravel()[mask_to_sq_rho],  dims = [[2]*(n), [2]*(n)])
    return(rho)




def prepare_from_rho_vec_to_rho_sq(n):
    s1 = np.array([['0', '1'],
                   ['2', '3']])
    s = s1
    for i in range(0, n-1):
        s = tensor_string_arrays(s, s1)
        
    mask_to_sq_rho = from_string_to_num(s, 4)
    T = tensor([Qobj(T_matrix)] * n)
    return(T, mask_to_sq_rho)
    

def tensor_string_arrays(s1, s2):
    n1 = len(s1[0])
    m1 = len(s1)
    n2 = len(s2[0])
    m2 = len(s2)
    m = m1*m2
    n = n1*n2
    result = [[0] * n for i in range(m)]
    for i in range(0, m):
        for j in range(0, n):
            result[i][j] = s1[i//m2][j//n2] + s2[i%m2][j%n2]
    return result


def from_string_to_num(l, n):
    res = np.zeros((len(l), len(l[0])))
    for i in range(0, len(l)):
        for j in range(0, len(l[0])):
            for k in range(len(l[i][j])):
                res[i, j] += int(l[i][j][-1-k])*n**(k)
    res = mask_to_sq_rho_to_int(res)
    return res


def mask_to_sq_rho_to_int(mask_to_sq_rho):
    mask_to_sq_rho = map(lambda x: list(map(lambda y: int(y), x)), mask_to_sq_rho)
    mask_to_sq_rho = np.asarray(list(mask_to_sq_rho))
    return mask_to_sq_rho