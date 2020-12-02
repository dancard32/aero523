import numpy as np
from numpy import linalg as LA
from readgri import writegri
from edgehash import edgehash
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
    #ind = int(np.ceil(flags.shape[0] * 0.1))
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

def genUE(u, Vcopy, V, E, IE, BE):
    Ecopy = E.copy(); Ucopy = u.copy()
    for i in range(Ecopy.shape[0]):
        n1, n2, n3 = Ecopy[i,:]
        x1 = V[int(n1),:]; x2 = V[int(n2),:]; x3 = V[int(n3),:]
        vals = np.array([(x2-x1)/2 +x1, (x3-x1)/2 +x1, (x3-x2)/2 +x2])
        
        # Generate nodes for each element
        nodes = np.array([])
        for k in vals:
            check, ind = vert_ind(Vcopy, k)
            if check:
                nodes = np.append(nodes, ind)
        
        if nodes.shape[0] == 3:
            # Ensure that the nodes are CCW
            if isCCW(Vcopy[int(nodes[0]),:], Vcopy[int(nodes[1]),:], Vcopy[int(nodes[2]),:]) != 1:
                nodes = np.flip(nodes)
            
            nodeint = np.array([n1, n2, n3])
            if isCCW(Vcopy[int(nodeint[0]),:], Vcopy[int(nodeint[1]),:], Vcopy[int(nodeint[2]),:]) != 1:
                nodeint = np.flip(nodeint)

            # Loop through the nodes
            for k in range(3):
                # Start at the nodes N1 -> N2 for consistency
                if Vcopy[int(nodes[k]),0] == vals[0,0] and Vcopy[int(nodes[k]),1] == vals[0,1]:
                    Ecopy[i,:] = np.array([n1, nodes[k], nodes[(k+2)%3]])   # Replace the ith element with new element

                    ind1 = np.array([nodes[k], n2, nodes[(k+1)%3]])
                    ind2 = np.array([nodes[k], nodes[(k+1)%3], nodes[(k+2)%3]])
                    ind3 = np.array([nodes[(k+1)%3], nodes[(k+2)%3], n3])
                    if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
                        ind1 = np.flip(ind1)
                    if isCCW(Vcopy[int(ind2[0]),:], Vcopy[int(ind2[1]),:], Vcopy[int(ind2[2]),:]) != 1:
                        ind2 = np.flip(ind2)
                    if isCCW(Vcopy[int(ind3[0]),:], Vcopy[int(ind3[1]),:], Vcopy[int(ind3[2]),:]) != 1:
                        ind3 = np.flip(ind3)
                    # Append new elements
                    Ecopy = np.append(Ecopy, np.transpose(np.array([[ind1[0]], [ind1[1]], [ind1[2]]])), axis=0)
                    Ecopy = np.append(Ecopy, np.transpose(np.array([[ind2[0]], [ind2[1]], [ind2[2]]])), axis=0)
                    Ecopy = np.append(Ecopy, np.transpose(np.array([[ind3[0]], [ind3[1]], [ind3[2]]])), axis=0)

                    for l in range(3):
                        Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)
                    break
            
        elif nodes.shape[0] == 2:
            node_ind = np.array([n1, n2, n3])
            
            if isCCW(Vcopy[int(node_ind[0]),:], Vcopy[int(node_ind[1]),:], Vcopy[int(node_ind[2]),:]) != 1:
                node_ind = np.flip(node_ind)
            
            ccw_count = 0
            for k in range(3):
                if isCCW(Vcopy[int(node_ind[k]),:], Vcopy[int(nodes[0]),:], Vcopy[int(nodes[1]),:]) != 1:
                    ccw_count += 1
            if ccw_count == 2:
                nodes = np.flip(nodes)

            for k in range(3):
                xn1 = Vcopy[int(node_ind[k]),:]; xn2 = Vcopy[int(node_ind[(k+1)%3]),:]; xn3 = Vcopy[int(node_ind[(k-1)%3]),:]

                test1 = (xn2-xn1)/2 + xn1; test2 = (xn1-xn3)/2 + xn3

                if '%.5f'%test1[0] == '%.5f'%Vcopy[int(nodes[1]),0] and '%.5f'%test1[1] == '%.5f'%Vcopy[int(nodes[1]),1] and '%.5f'%test2[0] == '%.5f'%Vcopy[int(nodes[0]),0] and '%.5f'%test2[1] == '%.5f'%Vcopy[int(nodes[0]),1]:
                    ind1 = np.array([node_ind[k], nodes[1], nodes[0]])
                    newnodes = np.array([node_ind[(k+1)%3], node_ind[(k+2)%3]])
            
            
            vn1 = Vcopy[int(nodes[0]),:]; vn2 = Vcopy[int(nodes[1]),:]
            xn1 = Vcopy[int(newnodes[0]),:]; xn2 = Vcopy[int(newnodes[1]),:]

            temp1 = (vn1-xn1)/2 + xn1; temp2 = (vn2-xn1)/2 + xn1
            theta1 = np.arccos(np.dot(temp1, temp2)/(LA.norm(temp1)*LA.norm(temp2)))

            temp1 = (vn1-xn2)/2 + xn2; temp2 = (vn2-xn2)/2 + xn2
            theta2 = np.arccos(np.dot(temp1, temp2)/(LA.norm(temp1)*LA.norm(temp2)))
                
            if theta1 >= theta2:
                ind2 = np.array([newnodes[1], nodes[0], nodes[1]])
                ind3 = np.array([newnodes[0], newnodes[1], nodes[1]])
            else:
                ind2 = np.array([newnodes[0], nodes[0], nodes[1]])
                ind3 = np.array([newnodes[0], newnodes[1], nodes[0]])
            
                            
            if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
                ind1 = np.flip(ind1)
            if isCCW(Vcopy[int(ind2[0]),:], Vcopy[int(ind2[1]),:], Vcopy[int(ind2[2]),:]) != 1:
                ind2 = np.flip(ind2)
            if isCCW(Vcopy[int(ind3[0]),:], Vcopy[int(ind3[1]),:], Vcopy[int(ind3[2]),:]) != 1:
                ind3 = np.flip(ind3)

            Ecopy[i,:] = np.array([ind1[0], ind1[1], ind1[2]])

            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind2[0]], [ind2[1]], [ind2[2]]])), axis=0)
            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind3[0]], [ind3[1]], [ind3[2]]])), axis=0)

            for k in range(2):
                Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)

        elif nodes.shape[0] == 1:
            
            for k in range(3):
                if vals[k,0] == Vcopy[int(nodes[0]),0] and vals[k,1] == Vcopy[int(nodes[0]),1]:
                    if k == 0:
                        ind1 = np.array([n1, nodes[0], n3])
                        ind2 = np.array([n2, n3, nodes[0]])
                    elif k == 1:
                        ind1 = np.array([n1, nodes[0], n2])
                        ind2 = np.array([n3, n2, nodes[0]])
                    elif k == 2:
                        ind1 = np.array([n2, nodes[0], n1])
                        ind2 = np.array([n3, n1, nodes[0]])
                    
            if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
                ind1 = np.flip(ind1)
            if isCCW(Vcopy[int(ind2[0]),:], Vcopy[int(ind2[1]),:], Vcopy[int(ind2[2]),:]) != 1:
                ind2 = np.flip(ind2)

            Ecopy[i,:] = np.array([ind1[0], ind1[1], ind1[2]])
            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind2[0]], [ind2[1]], [ind2[2]]])), axis=0)

            Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)

    for i in range(Ecopy.shape[0]):
        n1, n2, n3 = Ecopy[i,:]
        ind1 = np.array([n1, n2, n3])
        if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
            Ecopy[i,:] = np.array([Ecopy[i,2], Ecopy[i,1], Ecopy[i,0]])


    Ecopy = Ecopy.astype(int)
    return Ucopy, Ecopy

def genB(u, V, Vcopy, BE):
    Bcopy = BE.copy()
    for i in range(Bcopy.shape[0]):
        n1, n2, e1, bgroup = BE[i,:]
        xl = V[n1,:]; xr = V[n2,:]
        check, ind = vert_ind(Vcopy, 0.5*(xr+xl))
        if check:
            Bcopy[i,:] = np.array([n1, ind.item(), i, bgroup])
            Bcopy = np.append(Bcopy, np.transpose(np.array([[ind.item()],[n2], [Bcopy.shape[0]+1], [bgroup]])), axis=0)

    B0 = np.array([[-1,-1]]); B1 = B0.copy(); B2 = B0.copy(); B3 = B0.copy()
    for i in range(Bcopy.shape[0]):
        n1, n2, e, bname = Bcopy[i,:]
        if bname == 0:
            B0 = np.append(B0, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 1:
            B1 = np.append(B1, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 2:
            B2 = np.append(B2, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 3:
            B3 = np.append(B3, np.transpose(np.array([[n2], [n1]])), axis=0)
    B0 = B0[1:,:]; B1 = B1[1:,:]; B2 = B2[1:,:];  B3 = B3[1:,:];  
    B = [B0.astype(int), B1.astype(int), B2.astype(int), B3.astype(int)]

    return B

def isCCW(a, b, c):
    cross_val = (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])
    
    if cross_val > 0:
        cross_val = 1
    elif cross_val < 0:
        cross_val = -1
    else:
        cross_val = 0

    return cross_val

def vert_ind(Vvec, x):
    check = False; ind = np.array([])
    # Loop over the vertices
    for i in range(Vvec.shape[0]):
        # If this vertex exists return False
        
        if x[0] == Vvec[i,0] and x[1] == Vvec[i,1]:
            check = True; ind = np.append(ind, [i])
    
    return check, ind

def correctIE(IE):
    for i in range(IE.shape[0]):
        IE[i,:] = np.array([IE[i,0]-2, IE[i,1]-2, IE[i,2], IE[i,3]])

    return IE

def adapt(u, mach, V, E, IE, BE, filepath):

    flags = genflags(u, mach, V, E, IE, BE)         # Flag edges along the interior and exterior
    Vcopy = genV(flags, V, E, IE, BE)               # With the flags determine the nodes on the elements to split
    Ucopy, Ecopy = genUE(u, Vcopy, V, E, IE, BE)    # Determine the new Elements and with them the new U
    B = genB(u, V, Vcopy, BE)                       # With the old boundary edges determine which are on the edges
    IEcopy, BEcopy = edgehash(Ecopy, B)
    
    # Prepare for input to writegri       
    Mesh = {'V':Vcopy, 'E':Ecopy, 'IE':IEcopy, 'BE':BEcopy, 'Bname':['Engine', 'Exit', 'Outflow', 'Inflow'] }
    writegri(Mesh, filepath)

    return Ucopy, Vcopy, Ecopy, IEcopy, BEcopy
    
def plotmesh(V, BE, E):

    f = plt.figure(figsize=(10,10))
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    plt.scatter(V[:,0], V[:,1])
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='blue')
    plt.axis('equal'); plt.axis('off')
    f.tight_layout(); 
    plt.show()
    