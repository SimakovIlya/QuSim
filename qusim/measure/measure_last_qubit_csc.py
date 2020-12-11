from .projectors import *
from random import *


# R_NOT = tensor([Qobj(np.identity(4))] * 9 + [NOT_gate()])
# R = Proj0_gate(n = 10, p = 1, cur = 9)
# Proj0_gate_last_qubit_csc = csc_matrix(np.real(R.data).tocsc(), dtype = float)
# R = R_NOT * Proj1_gate(n = 10, p = 1, cur = 9)
# R_NOT_Proj1_gate_last_qubit_csc = csc_matrix(np.real(R.data).tocsc(), dtype = float)

Proj0_gate_last_qubit_csc = [0] * 10
R_NOT_Proj1_gate_last_qubit_csc = [0] * 10
for i in range(10):
    R_NOT = tensor([Qobj(np.identity(4))] * i + [NOT_gate()])
    R = Proj0_gate(n = i+1, p = 1, cur = i)
    Proj0_gate_last_qubit_csc[i] = csc_matrix(np.real(R.data).tocsc(), dtype = float)
    R = R_NOT * Proj1_gate(n = i+1, p = 1, cur = i)
    R_NOT_Proj1_gate_last_qubit_csc[i] = csc_matrix(np.real(R.data).tocsc(), dtype = float)


def measure_last_qubit_csc(rhoP_csc, eta = np.array([[1,0,0],[0,0,0]])): # input is rho in Pauli basis, return last qubit in |0>
    '''
    Measures the last (9, counting from 0) qubit and returns it in |0> state
    
    Args:
        rhoP (Qobj): rho-vector_csc in Pauli basis
        cur (int): qubit to measure
        
    Returns:
        rhoP (Qobj): rho-vector_csc in Pauli basis, last qubit in |0>
        m (int): 1 or -1, result of the measurement
    '''
    
    n = int(np.log2(rhoP_csc.shape[0])/2)
    # if (n != 10):
    #     print('measure_last_qubit: rewrite me!')

    P0 = (2**(n)*rhoP_csc[3, 0] + 1)/2
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
            rhoP_csc = (R_NOT_Proj1_gate_last_qubit_csc[n-1] * rhoP_csc) / P1
        else:
            rhoP_csc = (Proj0_gate_last_qubit_csc[n-1] * rhoP_csc) / P0
    elif p_out <= eta[mi, 1]:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = 1
        rhoP_csc = (R_NOT_Proj1_gate_last_qubit_csc[n-1] * rhoP_csc) / P1
    elif p_out <= eta[mi, 2]:
        m = -1
        if P0 == 0:
            print('measure_last_qubit P0')
            rhoP_csc = (R_NOT_Proj1_gate_last_qubit_csc[n-1] * rhoP_csc) / P1
        else:
            rhoP_csc = (Proj0_gate_last_qubit_csc[n-1] * rhoP_csc) / P0
    else:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = -1
        rhoP_csc = (R_NOT_Proj1_gate_last_qubit_csc[n-1] * rhoP_csc) / P1

    return(rhoP_csc, m)