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
        if '%.5f'%x[0] == '%.5f'%Vvec[i,0] and '%.5f'%x[1] == '%.5f'%Vvec[i,1]:
            check = False
            break
    
    return check

def isCCW(a, b, c):
    # Princeton implementation of CCW logical
    cross_val = (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])
    
    # Conditional statements
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
        
        if '%.5f'%x[0] == '%.5f'%Vvec[i,0] and '%.5f'%x[1] == '%.5f'%Vvec[i,1]:
            check = True; ind = np.append(ind, np.array([i]))
    
    return check, ind

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
        # Grab the node values of each given element
        n1, n2, n3 = Ecopy[i,:]
        x1 = V[int(n1),:]; x2 = V[int(n2),:]; x3 = V[int(n3),:]
        vals = np.array([(x2-x1)/2 +x1, (x3-x1)/2 +x1, (x3-x2)/2 +x2])
        
        # Generate nodes for each element
        nodes = np.array([])
        for k in vals:
            check, ind = vert_ind(Vcopy, k)
            if check:
                nodes = np.append(nodes, ind)
        
        # If three flags have been flagged
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
                if '%.5f'%Vcopy[int(nodes[k]),0] == '%.5f'%vals[0,0] and '%.5f'%Vcopy[int(nodes[k]),1] == '%.5f'%vals[0,1]:
                    Ecopy[i,:] = np.array([n1, nodes[k], nodes[(k+2)%3]])   # Replace the ith element with new element

                    # Start the indices based from generalized case
                    ind1 = np.array([nodes[k], n2, nodes[(k+1)%3]])
                    ind2 = np.array([nodes[k], nodes[(k+1)%3], nodes[(k+2)%3]])
                    ind3 = np.array([nodes[(k+1)%3], nodes[(k+2)%3], n3])

                    # Ensure that the nodes are CCW
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

                    # Append values of U to the updated U for initial starting condition
                    for l in range(3):
                        Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)
                    break
            
        # If two edges have been flagged
        elif nodes.shape[0] == 2:
            node_ind = np.array([n1, n2, n3])
            
            # Ensure CCW order
            if isCCW(Vcopy[int(node_ind[0]),:], Vcopy[int(node_ind[1]),:], Vcopy[int(node_ind[2]),:]) != 1:
                node_ind = np.flip(node_ind)
            
            # Ensure CCW order
            ccw_count = 0
            for k in range(3):
                if isCCW(Vcopy[int(node_ind[k]),:], Vcopy[int(nodes[0]),:], Vcopy[int(nodes[1]),:]) != 1:
                    ccw_count += 1
            if ccw_count == 2:
                nodes = np.flip(nodes)

            # Determine which node can be omitted
            for k in range(3):
                xn1 = Vcopy[int(node_ind[k]),:]; xn2 = Vcopy[int(node_ind[(k+1)%3]),:]; xn3 = Vcopy[int(node_ind[(k-1)%3]),:]
                test1 = (xn2-xn1)/2 + xn1; test2 = (xn1-xn3)/2 + xn3

                # Conditional to determine the starting nodes for triangle
                if '%.5f'%test1[0] == '%.5f'%Vcopy[int(nodes[1]),0] and '%.5f'%test1[1] == '%.5f'%Vcopy[int(nodes[1]),1] and '%.5f'%test2[0] == '%.5f'%Vcopy[int(nodes[0]),0] and '%.5f'%test2[1] == '%.5f'%Vcopy[int(nodes[0]),1]:
                    ind1 = np.array([node_ind[k], nodes[1], nodes[0]])
                    newnodes = np.array([node_ind[(k+1)%3], node_ind[(k+2)%3]])
            
            # Initialize values for test
            vn1 = Vcopy[int(nodes[0]),:]; vn2 = Vcopy[int(nodes[1]),:]
            xn1 = Vcopy[int(newnodes[0]),:]; xn2 = Vcopy[int(newnodes[1]),:]

            # Determine the angle of the first corner
            temp1 = (vn1-xn1)/2 + xn1; temp2 = (vn2-xn1)/2 + xn1
            theta1 = np.arccos(np.dot(temp1, temp2)/(LA.norm(temp1)*LA.norm(temp2)))

            # Determine the angle of the second corner
            temp1 = (vn1-xn2)/2 + xn2; temp2 = (vn2-xn2)/2 + xn2
            theta2 = np.arccos(np.dot(temp1, temp2)/(LA.norm(temp1)*LA.norm(temp2)))
                
            # Conditional to ensure no edge overlaps
            if theta1 >= theta2:
                # Re-arrangement of indices to ensure connectivity
                ind2 = np.array([newnodes[1], nodes[0], nodes[1]])
                ind3 = np.array([newnodes[0], newnodes[1], nodes[1]])
            else:
                # Re-arrangement of indices to ensure connectivity
                ind2 = np.array([newnodes[0], nodes[0], nodes[1]])
                ind3 = np.array([newnodes[0], newnodes[1], nodes[0]])
            
            # Ensure that they are CCW
            if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
                ind1 = np.flip(ind1)
            if isCCW(Vcopy[int(ind2[0]),:], Vcopy[int(ind2[1]),:], Vcopy[int(ind2[2]),:]) != 1:
                ind2 = np.flip(ind2)
            if isCCW(Vcopy[int(ind3[0]),:], Vcopy[int(ind3[1]),:], Vcopy[int(ind3[2]),:]) != 1:
                ind3 = np.flip(ind3)

            # Overwrite E, and append to E
            Ecopy[i,:] = np.array([ind1[0], ind1[1], ind1[2]])

            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind2[0]], [ind2[1]], [ind2[2]]])), axis=0)
            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind3[0]], [ind3[1]], [ind3[2]]])), axis=0)

            # Append to U for initial starting guess
            for k in range(2):
                Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)

        # If ones edges has been flagged
        elif nodes.shape[0] == 1:
            # Determine the starting node
            for k in range(3):
                if '%.5f'%vals[k,0] == '%=.5f'%Vcopy[int(nodes[0]),0] and '%.5f'%vals[k,1] == '%.5f'%Vcopy[int(nodes[0]),1]:
                    # Rearrange the nodes according to the node value
                    if k == 0:
                        ind1 = np.array([n1, nodes[0], n3])
                        ind2 = np.array([n2, n3, nodes[0]])
                    elif k == 1:
                        ind1 = np.array([n1, nodes[0], n2])
                        ind2 = np.array([n3, n2, nodes[0]])
                    elif k == 2:
                        ind1 = np.array([n2, nodes[0], n1])
                        ind2 = np.array([n3, n1, nodes[0]])
                    
            # Ensure CCW orientation
            if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
                ind1 = np.flip(ind1)
            if isCCW(Vcopy[int(ind2[0]),:], Vcopy[int(ind2[1]),:], Vcopy[int(ind2[2]),:]) != 1:
                ind2 = np.flip(ind2)

            # Overwrite and append E
            Ecopy[i,:] = np.array([ind1[0], ind1[1], ind1[2]])
            Ecopy = np.append(Ecopy, np.transpose(np.array([[ind2[0]], [ind2[1]], [ind2[2]]])), axis=0)

            # Append to U for initial start
            Ucopy = np.append(Ucopy, np.transpose(np.array([[u[i,0]], [u[i,1]], [u[i,2]], [u[i,3]]])), axis=0)

    # Double check to make sure CCW orientation
    for i in range(Ecopy.shape[0]):
        n1, n2, n3 = Ecopy[i,:]
        ind1 = np.array([n1, n2, n3])
        if isCCW(Vcopy[int(ind1[0]),:], Vcopy[int(ind1[1]),:], Vcopy[int(ind1[2]),:]) != 1:
            Ecopy[i,:] = np.array([Ecopy[i,2], Ecopy[i,1], Ecopy[i,0]])

    # Return E (as an integer array)
    Ecopy = Ecopy.astype(int)

    return Ucopy, Ecopy

def genB(u, V, Vcopy, BE):
    Bcopy = BE.copy()
    for i in range(Bcopy.shape[0]):
        # Node locations of boundary edges
        n1, n2, e1, bgroup = BE[i,:]
        xl = V[n1,:]; xr = V[n2,:]

        # Call function to determine if the vertex exists
        check, ind = vert_ind(Vcopy, (xl+xr)/2)
        if check:
            # If this vertex exists re-write B
            Bcopy[i,:] = np.array([n1, ind[0], e1, bgroup])
            Bcopy = np.append(Bcopy, np.transpose(np.array([[ind[0]],[n2], [Bcopy.shape[0]+1], [bgroup]])), axis=0)

    # Re-arrange B for input to edgehash()
    B0 = np.array([[-1,-1]]); B1 = B0.copy(); B2 = B0.copy(); B3 = B0.copy()
    for i in range(Bcopy.shape[0]):
        # Node vales
        n1, n2, e, bname = Bcopy[i,:]
        # Given the values of Bname append to the corresponding group
        if bname == 0:
            B0 = np.append(B0, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 1:
            B1 = np.append(B1, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 2:
            B2 = np.append(B2, np.transpose(np.array([[n1], [n2]])), axis=0)
        if bname == 3:
            B3 = np.append(B3, np.transpose(np.array([[n1], [n2]])), axis=0)
    
    # Output B#'s to B for input to edgehash()
    B0 = B0[1:,:]; B1 = B1[1:,:]; B2 = B2[1:,:];  B3 = B3[1:,:];  
    B = [B0.astype(int), B1.astype(int), B2.astype(int), B3.astype(int)]

    return B

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

    