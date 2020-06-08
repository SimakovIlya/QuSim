import sys
for i in range(len(sys.path)):
    if sys.path[i][len(sys.path[i])-len('qusim'):len(sys.path[i])] == 'qusim':
        sys.path.append(sys.path[i] + '/Pauli_basis')
        sys.path.append(sys.path[i] + '/gates')
        sys.path.append(sys.path[i] + '/measure')
        sys.path.append(sys.path[i] + '/plotting')


from import_Pauli_basis import *
from import_gates import *
from import_measure import *
from import_plotting import *