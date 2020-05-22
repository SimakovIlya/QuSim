from projectors import *
from random import *

N_tmp = 10
Proj0_gate_for_measure_qubit_csc_func = [csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 0).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 1).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 2).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 3).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 4).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 5).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 6).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 7).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 8).data).tocsc(), dtype = float),
                                         csc_matrix(np.real(Proj0_gate(n = N_tmp, p = 1, cur = 9).data).tocsc(), dtype = float)]
qc_tmp = QCircuit(N_tmp)
not_proj1_tmp = to_Pauli_T_matrix(np.array([[0, 1], [1, 0]]) @ np.array([[0, 0], [0, 1]]))
R_NOT_Proj1_gate_for_measure_qubit_csc_func = [csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 0)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 1)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 2)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 3)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 4)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 5)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 6)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 7)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 8)).data).tocsc(), dtype = float),
                                               csc_matrix(np.real((qc_tmp.make_gate(not_proj1_tmp, 9)).data).tocsc(), dtype = float)]



def measure_qubit_csc(rhoP_csc, cur, eta = np.array([[1,0,0],[0,0,0]])):
    '''
    Measures qubit and returns it in |0> state
    
    Args:
        rhoP_csc: rho-vector_csc in Pauli basis
        cur (int): qubit to measure
        eta (ndarray): erros)
        
    Returns:
        rhoP_csc: rho-vector_csc in Pauli basis, last qubit in |0>
        m (int): 1 or -1, result of the measurement
    '''
    n = int(np.log2(rhoP_csc.shape[0])/2)
    P0 = (2**(n)*rhoP_csc[3 * 4**(n - cur - 1), 0] + 1)/2
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
        rhoP_csc = (Proj0_gate_for_measure_qubit_csc_func[cur] * rhoP_csc) / P0
    elif p_out <= eta[mi, 1]:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = -1
        rhoP_csc = (R_NOT_Proj1_gate_for_measure_qubit_csc_func[cur] * rhoP_csc) / P1
    elif p_out <= eta[mi, 2]:
        m = 1
        if P0 == 0:
            print('measure_last_qubit P0')
        R = Proj0_gate(n = n, p = 1, cur = cur)
        rhoP_csc = (Proj0_gate_for_measure_qubit_csc_func[cur] * rhoP_csc) / P0
    else:
        if P1 == 0:
            print('measure_last_qubit P1')
        m = -1
        rhoP_csc = (R_NOT_Proj1_gate_for_measure_qubit_csc_func[cur] * rhoP_csc) / P1

    return(rhoP_csc, m)