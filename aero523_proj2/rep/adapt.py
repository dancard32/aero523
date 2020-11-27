import numpy as np
from numpy import linalg as LA
from readgri import edgehash, writegri

def mach_perp(u, nhat):
    uvel = u[1]/u[0]; v = u[2]/u[0]

    q = np.dot(np.array([uvel, v]), nhat)

    P = (1.4 - 1)*(u[3] - 0.5*u[0]*q**2)
    H = (u[3] + P)/u[0]
    c = np.sqrt(0.4*(H - 0.5*q**2))
    mach = q/c

    return mach

def adapt(u, mach, V, E, IE, BE):

    flags = np.zeros(((IE.shape[0] + BE.shape[0]),2)); flags[:,0] = np.arange(flags.shape[0]);k = 0
    for i in range(IE.shape[0]):
        n1, n2, e1, e2 = IE[i,:]
        xl = V[n1,:]; xr = V[n2,:]
        machl = mach[e1]; machr = mach[e2]
        dx = xr - xl; deltal = LA.norm(dx)
        eps = abs(machr - machl)*deltal
        
        flags[k,1] += eps; k += 1
    for i in range(BE.shape[0]):
        n1, n2, e1, bgroup = BE[i,:]
        #if bgroup == 0: # Engine
        xl = V[n1,:]; xr = V[n2,:]
        uedge = u[e1,:]
        dx = xr - xl; deltal = LA.norm(dx)
        nhat = np.array([dx[1], -dx[0]])/deltal
        machperp = mach_perp(uedge, nhat)
        eps = abs(machperp)*deltal
        
        flags[k,1] += eps
        k += 1

    # Sort from largest to smallest errors
    flags = flags[flags[:,1].argsort()]; flags = np.flipud(flags)
    # Remove all outliers to be refined
    ind = int(np.ceil(flags.shape[0] * 0.03))
    flags[ind:(flags.shape[0]-1),1] = 0
    flags = flags[flags[:,0].argsort()]

    print(flags)
    
    print(E)
    Vcopy = V.copy(); Ecopy = E.copy(); Bcopy = BE.copy(); Bcopy = Bcopy[:,0:2]; k = 0
    for i in range(IE.shape[0]):
        err = flags[k,1]
        if err > 0:
            n1, n2, e1, e2 = IE[i,:]
            
            n3, n4, n5 = E[e1,:]
            n6, n7, n8 = E[e2,:]

        k += 1
    for i in range(BE.shape[0]):
        err = flags[k,1]
        if err > 0:
            edge1, edge2, e1, ig = BE[i,:]
            print(e1, '\n',E[e1,:])
            
            n1, n2, n3 = E[e1,:]
            x1 = V[n1,:]; x2 = V[n2,:]; x3 = V[n3,:]

            # Conditionals to prevent duplicate nodes
            Vcopy = np.append(Vcopy, np.array([(x2-x1)/2 +x1]), axis=0)
            Vcopy = np.append(Vcopy, np.array([(x3-x1)/2 +x1]), axis=0)
            Vcopy = np.append(Vcopy, np.array([(x3-x2)/2 +x2]), axis=0)

            # Conditionals to determine the new boundaries


            # Create new elements

        k += 1

    print(Vcopy)