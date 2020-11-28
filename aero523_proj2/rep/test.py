import numpy as np
from numpy import linalg as LA
import random
import matplotlib.pyplot as plt
import time

# Project specific functions
from readgri import readgri, writegri
from plotmesh import plotmesh
from flux import RoeFlux
from fvm import solve
from adapt import adapt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def getIC(alpha, Ne):
    alpha = np.deg2rad(alpha); Minf = 2.2; gam = 1.4
    uinf = np.array([1, Minf*np.cos(alpha), Minf*np.sin(alpha), 1/(gam*(gam-1)) + Minf**2/2])

    u0 = np.zeros((Ne, 4))
    for i in range(4):
        u0[:,i] = uinf[i]
    u0[abs(u0) < 10**-10]

    return u0

def post_process(u):
    uvel = u[:,1]/u[:,0]; v = u[:,2]/u[:,0]
    
    q = np.zeros(u.shape[0])
    for i in range(u.shape[0]):
        q[i] = LA.norm(np.array([uvel[i], v[i]]))

    P = (1.4 - 1)*(u[:,3] - 0.5*u[:,0]*q**2)
    H = (u[:,3] + P)/u[:,0]
    c = np.sqrt(0.4*(H - 0.5*q**2))
    mach = q/c

    Pt = P*(1 + 0.5*0.4*mach**2)**(1.4/0.4)

    return mach, Pt

def _3by3_():
    mesh = readgri('test0.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    u = getIC(1, E.shape[0])
    test = np.random.rand(8)
    u[:,0] += test
    #u[:,0] += np.linspace(0,2, num=8)
    mach, Pt = post_process(u)

    adapt(u, mach, V, E, IE, BE)

def main():
    mesh = readgri('mesh0.gri')
    
    start = time.time()
    u, err, ATPR, V, E, BE, IE = solve(1, mesh); end = time.time(); print('Elapsed Time %.2f'%(end - start))
    mach, pt = post_process(u)

    adapt(u, mach, V, E, IE, BE)


if __name__ == "__main__":
    _3by3_()