import numpy as np
from numpy import linalg as LA
import random
import matplotlib.pyplot as plt
import time

# Project specific functions
from readgri import readgri, writegri
from flux import RoeFlux
from fvm import solve
from adapt import adapt

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

def FVMIC(alpha, Ne):
    alpha = np.deg2rad(alpha); Minf = 2.2; gam = 1.4
    uinf = np.array([1, Minf*np.cos(alpha), Minf*np.sin(alpha), 1/(gam*(gam-1)) + Minf**2/2])

    u0 = np.zeros((Ne, 4))
    for i in range(4):
        u0[:,i] = uinf[i]
    u0[abs(u0) < 10**-10]

    return u0

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

def plotmesh(V, BE, E):

    f = plt.figure(figsize=(10,10))
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    plt.scatter(V[:,0], V[:,1])
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='blue')
    plt.axis('equal'); plt.axis('off')
    f.tight_layout(); 
    plt.show()

def _3by3_():
    alpha = 1
    mesh = readgri('test0.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    
    u = getIC(alpha, E.shape[0])
    test = np.random.rand(8)
    test = np.array([0.53189414, 0.58794786, 0.4576574,  0.24801386, 0.72982654, 0.02918595, 0.19271085, 0.26531498]) # Double flag case
    test = np.array([0.0438649, 0.41647822, 0.98259851, 0.77928424, 0.53674341, 0.82151621, 0.0979844,  0.06814286])  # Double flag edge case
    test = np.array([0.69202655, 0.35823124, 0.49821811, 0.31857253, 0.58021842, 0.42569984, 0.94716888, 0.89860787]) # Double flag edge case   
    test = np.array([0.47834163, 0.99764052, 0.36555306, 0.24986927, 0.26026875, 0.34121071, 0.50115206, 0.04995461])  # Double flag edge case

    
    mesh = readgri('test0.gri'); E = mesh['E']
    u = FVMIC(alpha, E.shape[0])
    
    u, err, ATPR, V, E, BE, IE = solve(alpha, u, mesh)
    mach, pt = post_process(u)

    # Adapt the mesh
    u, V, E, IE, BE = adapt(u, mach, V, E, IE, BE, 'test' + str(0+1) + '.gri')

    plotmesh(V, BE, E)


def main():
    mesh = readgri('mesh0.gri'); E = mesh['E']
    u = FVMIC(1, E.shape[0])

    start = time.time()
    u, err, ATPR, V, E, BE, IE = solve(1, u, mesh); end = time.time(); print('Elapsed Time %.2f'%(end - start))
    mach, pt = post_process(u)

    u, V, E, IE, BE = adapt(u, mach, V, E, IE, BE, 1)


if __name__ == "__main__":
    _3by3_()
    #main()