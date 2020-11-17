import numpy as np
from readgri import readgri, writegri


def solve(adaptnum):
    alpha = 0; Minf = 2.2; gam = 1.4
    uinf = np.transpose(np.array([1, Minf*np.cos(alpha), Minf*np.sin(alpha), 1/(gam*(gam-1)) + Minf**2/2]))

    mesh = readgri('mesh' + adaptnum + '.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    u = np.zeros((E.shape[0], 4))
    

    return 0