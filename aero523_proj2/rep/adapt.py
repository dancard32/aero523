import numpy as np
from numpy import linalg as LA
from readgri import edgehash, writegri
import matplotlib.pyplot as plt

def mach_perp(u, nhat):
    uvel = u[1]/u[0]; v = u[2]/u[0]         # Calculate the velocity
    q = np.dot(np.array([uvel, v]), nhat)   # Determine the perpendicular speed
    P = (1.4 - 1)*(u[3] - 0.5*u[0]*q**2)    # Calculate pressure
    H = (u[3] + P)/u[0]                     # Calculate enthalpy
    c = np.sqrt(0.4*(H - 0.5*q**2))         # Calculate speed of sound
    mach = q/c                              # Calculate the Mach number 

    return mach

def check_vert(Vvec, x):
    check = True
    # Loop over the vertices
    for i in range(Vvec.shape[0]):
        # If this vertex exists return False
        if x[0] == Vvec[i,0] and x[1] == Vvec[i,1]:
            check = False
            break
    
    return check

def genflags(u, mach, V, E, IE, BE):

    # Pre-allocate flag array
    flags = np.zeros(((IE.shape[0] + BE.shape[0]),2)); flags[:,0] = np.arange(flags.shape[0]);k = 0

    # Iterate over the interior edges
    for i in range(IE.shape[0]):
        n1, n2, e1, e2 = IE[i,:]            # Nodes and elements from interior edge
        xl = V[n1,:]; xr = V[n2,:]          # Vertice values
        machl = mach[e1]; machr = mach[e2]  # Mach numbers at each element
        dx = xr - xl; deltal = LA.norm(dx)  # Determine the length of the edge
        eps = abs(machr - machl)*deltal     # Calculate the error
        
        flags[k,1] += eps; k += 1           # Add the edge error

    # Iterate over the boundary edges
    for i in range(BE.shape[0]):
        n1, n2, e1, bgroup = BE[i,:]        # Node and elements from boundary edge
        if bgroup == 0: # Engine
            xl = V[n1,:]; xr = V[n2,:]          # Vertice values
            uedge = u[e1,:]                     # State at edge
            dx = xr - xl; deltal = LA.norm(dx)  # Determine the length of the edge
            nhat = np.array([dx[1], -dx[0]])/deltal # Determine the normal off the boundary edge
            machperp = mach_perp(uedge, nhat)   # Calculate the perpendicular Mach number
            eps = abs(machperp)*deltal          # Calculate error
            
            flags[k,1] += eps                   # Add the edge error
        k += 1

    # Sort from largest to smallest errors
    flags = flags[flags[:,1].argsort()]; flags = np.flipud(flags)
    # Remove all outliers to be refined
    ind = int(np.ceil(flags.shape[0] * 0.03))
    #ind = int(np.ceil(flags.shape[0] * 0.65))
    flags[ind:(flags.shape[0]-1),1] = 0

    # Sort the errors increasing the edge number to iterate
    flags = flags[flags[:,0].argsort()]

    return flags

def genV(flags, V, E, IE, BE):
    Vcopy = V.copy(); k = 0
    for i in range(IE.shape[0]):
        err = flags[k,1]
        if err > 0:
            ig, ig, e1, e2 = IE[i,:]
            for j in np.array([e1,e2]):
                n1, n2, n3 = E[j,:]
                x1 = V[n1,:]; x2 = V[n2,:]; x3 = V[n3,:]
                
                # Conditionals to prevent duplicate nodes
                if check_vert(Vcopy, (x2-x1)/2 +x1):
                    Vcopy = np.append(Vcopy, np.array([(x2-x1)/2 +x1]), axis=0)
                if check_vert(Vcopy, (x3-x1)/2 +x1):
                    Vcopy = np.append(Vcopy, np.array([(x3-x1)/2 +x1]), axis=0)
                if check_vert(Vcopy, (x3-x2)/2 +x2):
                    Vcopy = np.append(Vcopy, np.array([(x3-x2)/2 +x2]), axis=0)
        k += 1
    
    for i in range(BE.shape[0]):
        err = flags[k,1]
        if err > 0:
            ig, ig, e1, ig = BE[i,:]
            n1, n2, n3 = E[e1,:]
            x1 = V[n1,:]; x2 = V[n2,:]; x3 = V[n3,:]

            # Conditionals to prevent duplicate nodes
            if check_vert(Vcopy, (x2-x1)/2 +x1):
                Vcopy = np.append(Vcopy, np.array([(x2-x1)/2 +x1]), axis=0)

            if check_vert(Vcopy, (x3-x1)/2 +x1):
                Vcopy = np.append(Vcopy, np.array([(x3-x1)/2 +x1]), axis=0)            

            if check_vert(Vcopy, (x3-x2)/2 +x2):
                Vcopy = np.append(Vcopy, np.array([(x3-x2)/2 +x2]), axis=0)
        k += 1

    return Vcopy

def isboundary(nodestate, BEvec, Vvec):
    check = False
    for k in range(3):
        node = Vvec[int(nodestate[k])]
        for i in range(BEvec.shape[0]):
            n1, ig, ig, ig = BEvec[i,:]        # Node and elements from boundary edge
            x1 = Vvec[n1,:]
            if node[0] == x1[0] and node[1] == x1[1]:
                check = True

    return check
def vert_ind(Vvec, x):
    check = False; ind = -1
    # Loop over the vertices
    for i in range(Vvec.shape[0]):
        # If this vertex exists return False
        if x[0] == Vvec[i,0] and x[1] == Vvec[i,1]:
            check = True; ind = i
            break
    
    return check, ind

def adapt(u, mach, V, E, IE, BE):

    flags = genflags(u, mach, V, E, IE, BE)
    Vcopy = genV(flags, V, E, IE, BE)

    Ecopy = E.copy()
    for i in range(Ecopy.shape[0]):
        n1, n2, n3 = Ecopy[i,:]
        x1 = V[int(n1),:]; x2 = V[int(n2),:]; x3 = V[int(n3),:]
        vals = np.array([(x2-x1)/2 +x1, (x3-x1)/2 +x1, (x3-x2)/2 +x2])

        nodes = np.array([])
        for k in vals:
            check, ind = vert_ind(Vcopy, k)
            if check:
                nodes = np.append(nodes, ind)

        if nodes.shape[0] > 2:
            
            Ecopy[i,:] = np.array([n1, nodes[1], nodes[0]])

            Ecopy = np.append(Ecopy, np.transpose(np.array([[nodes[2]], [n2], [nodes[0]]])), axis=0)
            Ecopy = np.append(Ecopy, np.transpose(np.array([[nodes[2]], [nodes[1]], [nodes[0]]])), axis=0)
            Ecopy = np.append(Ecopy, np.transpose(np.array([[nodes[1]], [nodes[2]], [n3]])), axis=0)


        elif nodes.shape[0] > 1:
            pass

        elif nodes.shape[0] > 0:            
            nodesort = np.array([n1,n2,n3])
            nodesort = nodesort[nodesort.argsort()]
            
            # Fix if the node rests at a boundary
            if isboundary(nodesort, BE, V):    
                Ecopy = np.append(Ecopy, np.transpose(np.array([[nodes[0]], [nodesort[0]], [nodesort[2]]])), axis=0)

                
            
    print(Vcopy.shape[0])
    print(Ecopy.shape[0])
    plotmesh(Vcopy, BE, Ecopy)
    
def plotmesh(V, B, E):

    f = plt.figure(figsize=(12,12))
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    plt.scatter(V[:,0], V[:,1])
    #for i in range(BE.shape[0]):
    #    plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='black')
    plt.axis('equal'); plt.axis('off')
    f.tight_layout(); 
    plt.show()