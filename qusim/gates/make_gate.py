from import_Pauli_basis import *


class QCircuit:
    def __init__(self, n):
        self.n = n


    def make_gate(self, T, cur, cur1 = -1):
        '''
        Returns trasfer matrix of the gate in n-dimensional space
        
        Args:
            T (Qobj): trasfer matrix of the gate in 1-dim space
            cur (int): cursor to the qubit on which the gate is applied
            cur1 (int): cursor to the second qubit on which the gate is applied
            
        Returns:
            Qobj: trasfer matrix (4^n x 4^n)
        '''
        assert(max(cur, cur1) <= self.n)
        if T.data.shape[0] == 4 and cur1 == -1:
            R = tensor([Qobj(np.identity(4))] * cur + [T] + [Qobj(np.identity(4))] * (self.n - 1 - cur))
        elif T.data.shape[0] == 16 and cur1 != -1:
            l = np.arange(2, self.n).tolist()
            l.insert(min(cur, cur1), 0)
            l.insert(max(cur1, cur), 1)
            R = tensor([T] + [Qobj(np.identity(4))] * (self.n - 2)).permute(l)
        else:
            print('make_gate ERROR')
        return(R)