from ..Pauli_basis import *
from ..gates import *
from random import *



def Proj0_gate(n, cur, p): 
    return tensor([Qobj(np.identity(4))] * cur + [Qobj(to_Pauli_T_matrix(Proj0(np.sqrt(p))))]\
                + [Qobj(np.identity(4))] * (n - 1 - cur))
    
def Proj1_gate(n, cur, p): 
    return tensor([Qobj(np.identity(4))] * cur + [Qobj(to_Pauli_T_matrix(Proj1(np.sqrt(p))))]\
                + [Qobj(np.identity(4))] * (n - 1 - cur))
    
def Projplus_gate(n, cur, p): 
    return tensor([Qobj(np.identity(4))] * cur + [Qobj(to_Pauli_T_matrix(Projplus(np.sqrt(p))))]\
                + [Qobj(np.identity(4))] * (n - 1 - cur))
    
def Projminus_gate(n, cur, p): 
    return tensor([Qobj(np.identity(4))] * cur + [Qobj(to_Pauli_T_matrix(Projminus(np.sqrt(p))))]\
                + [Qobj(np.identity(4))] * (n - 1 - cur))




def Proj0(p): # NOT in Pauli
    return(np.array([[1/p, 0],
                     [0, 0]]))
def Proj1(p): # NOT in Pauli
    return(np.array([[0, 0],
                     [0, 1/p]]))
def Projplus(p): # NOT in Pauli
    return(1/2/p * np.array([[1, 1],
                             [1, 1]]))
def Projminus(p): # NOT in Pauli
    return(1/2/p * np.array([[1, -1],
                            [-1, 1]]))