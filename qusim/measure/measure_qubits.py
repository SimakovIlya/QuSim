from projectors import *



Proj0_gate_for_all_data_qubits_func = [Proj0_gate(n = 1, p = 1, cur = 0),
                                       Proj0_gate(n = 2, p = 1, cur = 1),
                                       Proj0_gate(n = 3, p = 1, cur = 2),
                                       Proj0_gate(n = 4, p = 1, cur = 3),
                                       Proj0_gate(n = 5, p = 1, cur = 4),
                                       Proj0_gate(n = 6, p = 1, cur = 5),
                                       Proj0_gate(n = 7, p = 1, cur = 6),
                                       Proj0_gate(n = 8, p = 1, cur = 7),
                                       Proj0_gate(n = 9, p = 1, cur = 8)]
Proj1_gate_for_all_data_qubits_func = [Proj1_gate(n = 1, p = 1, cur = 0),
                                       Proj1_gate(n = 2, p = 1, cur = 1),
                                       Proj1_gate(n = 3, p = 1, cur = 2),
                                       Proj1_gate(n = 4, p = 1, cur = 3),
                                       Proj1_gate(n = 5, p = 1, cur = 4),
                                       Proj1_gate(n = 6, p = 1, cur = 5),
                                       Proj1_gate(n = 7, p = 1, cur = 6),
                                       Proj1_gate(n = 8, p = 1, cur = 7),
                                       Proj1_gate(n = 9, p = 1, cur = 8)]




def measure_all_data_qubits(rhoP, eta = np.array([[1,0,0],[0,0,0]])):
    '''
    Consistently measures all data qubit
    
    Args:
        rhoP (Qobj): rho-vector in Pauli basis
        
    Returns:
        M_data (ndarray): array of 1 or -1, results of the measurements
    '''
    n = int(np.log2(rhoP.shape[0])/2)
    # if (n != 9):
    #     print('func measure_all_data_qubits askes you think again and make 9 data qubits')
    M_data = np.zeros((n))
    for i in range(n, 0, -1):
        P0 = (2**(i)*rhoP[3, 0] + 1)/2
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
            rhoP = (Proj0_gate_for_all_data_qubits_func[i-1] * rhoP) / P0
        elif p_out <= eta[mi, 1]:
            if P1 == 0:
                print('measure_last_qubit P1')
            m = -1
            rhoP = (Proj1_gate_for_all_data_qubits_func[i-1] * rhoP) / P1
        elif p_out <= eta[mi, 2]:
            m = 1
            if P0 == 0:
                print('measure_last_qubit P0')
            rhoP = (Proj0_gate_for_all_data_qubits_func[i-1] * rhoP) / P0
        else:
            if P1 == 0:
                print('measure_last_qubit P1')
            m = -1
            rhoP = (Proj1_gate_for_all_data_qubits_func[i-1] * rhoP) / P1
        M_data[i - 1] = m
        if i > 1:
            rhoP = Qobj(2*rhoP[::4], dims = [[4]*(i-1), [1]*(i-1)])
    return(M_data) 