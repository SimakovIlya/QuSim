from projectors import *



# R_NOT = tensor([Qobj(np.identity(4))] * 9 + [NOT_gate()])
# Proj0_gate_last_qubit = Proj0_gate(n = 10, p = 1, cur = 9)
# R_NOT_Proj1_gate_last_qubit = R_NOT * Proj1_gate(n = 10, p = 1, cur = 9)

def measure_last_qubit(rhoP, eta = np.array([[1,0,0],[0,0,0]])): # input is rho in Pauli basis, return last qubit in |0>
    '''
    Measures the last (9, counting from 0) qubit and returns it in |0> state
    
    Args:
        rhoP (Qobj): rho-vector in Pauli basis
        ideal (bool): whether gate is ideal or not
        
    Returns:
        rhoP (Qobj): rho-vector in Pauli basis, last qubit in |0>
        m (int): 1 or -1, result of the measurement
    '''
    
    n = int(np.log2(rhoP.shape[0])/2)
    R_NOT = tensor([Qobj(np.identity(4))] * (n-1) + [NOT_gate()])
    Proj0_gate_last_qubit = Proj0_gate(n = n, p = 1, cur = n-1)
    R_NOT_Proj1_gate_last_qubit = R_NOT * Proj1_gate(n = n, p = 1, cur = n-1)
    # if (n != 10):
    #     print('measure_last_qubit: rewrite me!')
    P0 = (2**(n)*rhoP[3, 0] + 1)/2
    P1 = 1 - P0
    if random() <= P0:
        mi = 0
    else:
        mi = 1
    p_out = random()
    if p_out <= eta[mi, 0]:
        m = 1
        if P0 == 0:
            print('measure_last_qubit P0')
        rhoP = (Proj0_gate_last_qubit * rhoP) / P0
    elif p_out <= eta[mi, 1]:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = -1
        rhoP = (R_NOT_Proj1_gate_last_qubit * rhoP) / P1
    elif p_out <= eta[mi, 2]:
        m = 1
        if P0 == 0:
            print('measure_last_qubit P0')
        rhoP = (Proj0_gate_last_qubit * rhoP) / P0
    else:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = -1
        rhoP = (R_NOT_Proj1_gate_last_qubit * rhoP) / P1
        
    return(rhoP, m)