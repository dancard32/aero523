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

    flags = np.zeros((E.shape[0],3)); flags[:,0] = np.arange(E.shape[0])
    for i in range(IE.shape[0]):
        n1, n2, e1, e2 = IE[i,:]
        xl = V[n1,:]; xr = V[n2,:]
        machl = mach[e1]; machr = mach[e2]
        dx = xr - xl; deltal = LA.norm(dx)
        eps = abs(machr - machl)*deltal
        
        flags[e1,1] += 1; flags[e1,2] += eps
        flags[e2,1] += 1; flags[e2,2] += eps
    for i in range(BE.shape[0]):
        n1, n2, e1, bgroup = BE[i,:]
        #if bgroup == 0: # Engine
        xl = V[n1,:]; xr = V[n2,:]
        uedge = u[e1,:]
        dx = xr - xl; deltal = LA.norm(dx)
        nhat = np.array([dx[1], -dx[0]])/deltal
        machperp = mach_perp(uedge, nhat)
        eps = abs(machperp)*deltal
        
        flags[e1,1] += 1; flags[e1,2] += eps

    # Sort from largest to smallest errors
    flags = flags[flags[:,2].argsort()]; flags = np.flipud(flags)
    # Remove all outliers to be refined
    ind = int(np.ceil(E.shape[0] * 0.3))
    flags[ind:(E.shape[0]-1),1] = 0
    flags = flags[flags[:,0].argsort()]
    
    # Flag all edges - test
    flags[:,1] = 3
    
    Etemp = E.copy(); Vtemp = V.copy()
    for i in range(1):
        n1, n2, n3 = E[i,:]
        x1 = V[n1,:]; x2 = V[n2,:]; x3 = V[n3,:]


        print(n1,n2,n3)
        if flags[i,1] == 3:
            [n4,n5,n6] = np.arange(1,4) + i

            # Refine Elements
            Etemp[int(flags[i,0]),:] = np.array([n1,n5,n4])
            Etemp = np.insert(Etemp, n4, [[n5,n2,n6]], axis=0)
            Etemp = np.insert(Etemp, n5, [[n6,n3,n4]], axis=0)
            Etemp = np.insert(Etemp, n6, [[n5,n6,n4]], axis=0)

            # Refine Vertices
            Vtemp = np.insert(Vtemp, [n4,n5,n6], [x3 - x1], axis=0)
            Vtemp = np.insert(Vtemp, [n4,n5,n6], [x2 - x1], axis=0)
            Vtemp = np.insert(Vtemp, [n4,n5,n6], [x3 - x2], axis=0)
    print(np.shape(Etemp))
    print(np.shape(Vtemp))

    


    B = BE[:,0:2]; 
    IE, BE = edgehash(E,B)





    """
    # Loop over values that need refinement
    for i in range(flags.shape[0]):
        if flags[i,2] == 0: # No refinement needed
            break
        else:
            if flags[i,3] == 10:
                [n1,n2,n3] = E[int(flags[i,0]),:]
                x1 = V[n1,:]; x2 = V[n2,:]; x3 = V[n3,:]

                if flags[i,1] == 3:
                    [n4,n5,n6] = np.arange(1,4) + E.shape[0]

                    # Refine Elements
                    E[int(flags[i,0]),:] = np.array([n1,n5,n4])
                    E = np.append(E, [[n5,n2,n6]], axis=0)
                    E = np.append(E, [[n6,n3,n4]], axis=0)
                    E = np.append(E, [[n5,n6,n4]], axis=0)

                    # Refine Vertices
                    V = np.append(V, [x3 - x1], axis=0)
                    V = np.append(V, [x2 - x1], axis=0)
                    V = np.append(V, [x3 - x2], axis=0)

                elif flags[i,1] == 2:
                    [n4,n5] = np.arange(1,3) + E.shape[0]

                    # Refine Elements
                    E[int(flags[i,0]),:] = np.array([n1,n2,n4])
                    E = np.append(E, [[n2,n5,n4]], axis=0)
                    E = np.append(E, [[n5,n3,n4]], axis=0)

                    # Refine Vertices
                    V = np.append(V, [x3 - x1], axis=0)
                    V = np.append(V, [x3 - x2], axis=0)
                elif flags[i,1] == 1:
                    n4 = E.shape[0] + 1

                    # Refine Elements
                    E[int(flags[i,0]),:] = np.array([n1,n4,n3])
                    E = np.append(E, [[n4,n2,n3]], axis=0)

                    # Refine Vertices
                    V = np.append(V, [x2 - x1], axis=0)
    
    """