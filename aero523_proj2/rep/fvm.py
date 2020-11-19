import numpy as np
from numpy import linalg as LA
from readgri import readgri, writegri

from flux import RoeFlux

def getIC(alpha, Ne):
    alpha = np.deg2rad(alpha); Minf = 2.2; gam = 1.4
    uinf = np.array([1, Minf*np.cos(alpha), Minf*np.sin(alpha), 1/(gam*(gam-1)) + Minf**2/2])

    u0 = np.zeros((Ne, 4))
    for i in range(Ne):
        u0[i,:] = uinf

    return u0

def solve():
    mesh = readgri('mesh0.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    # Get the initial state
    u0 = getIC(0, E.shape[0]); u = u0.copy()

    # Initialize the residual matrix
    R = np.zeros((E.shape[0], 4)); dta = R.copy(); err = np.array([0])
    
    print('Residual Norm\n' + 50*'-'); 
    # Iterate while out of tolerance
    for itr in range(5):
        R *= 0; dta *= 0    # Re-initialize

        # Loop over Interior Edges
        for i in range(IE.shape[0]):
            n1, n2, e1, e2 = IE[i,:]
            x1 = V[n1,:]; x2 = V[n2,:]
            
            u1 = u[e1,:]; u2 = u[e2,:]
            dx = x2 - x1; deltal = LA.norm(dx); 
            nhat = np.array([x2[1] - x1[1], x2[0] - x1[0]])/deltal

            # Determine the flux
            F, ls = RoeFlux(u1, u2, nhat, False)
            R[e1,:] -= F*deltal; R[e2,:] += F*deltal
            dta[e1,:] += abs(max(ls))*deltal; dta[e2,:] += abs(max(ls))*deltal

        # Implement Boundaries
        for i in range(BE.shape[0]):
            n1, n2, elem, bgroup = BE[i,:]
            x1 = V[n1,:]; x2 = V[n2,:]
            
            uedge = u[elem,:]
            dx = x2 - x1; deltal = LA.norm(dx); 
            nhat = np.array([x2[1] - x1[1], x2[0] - x1[0]])/deltal

            if bgroup == 0:                                 # Engine

                # Inviscid Boundary Conditions
                vplus = np.array([uedge[1], uedge[2]])/uedge[0]
                vb = vplus - np.dot(vplus, nhat)*nhat
                rhob = (1.4-1)*(uedge[3] - 0.5*uedge[0]*LA.norm(vb))

                F, ls = RoeFlux(uedge, uedge, nhat, False) 
                F = np.array([0, rhob*nhat[0], rhob*nhat[1], 0])
                R[elem,:] -= F*deltal

            elif bgroup == 2 or bgroup ==1:                 # Outflow

                # Down-stream conditions are equal
                F, ls = RoeFlux(uedge, uedge, nhat, False) 
                R[elem,:] -= F*deltal

            elif bgroup == 3:                               # Inflow
                
                # Set to the initial state for inflow
                F, ls = RoeFlux(uedge, u0[0,:], nhat, False) 
                R[elem,:] += F*deltal                    
            
            dta[elem,:] += abs(max(ls))*deltal

        dta *= 2/dta
        u -= np.multiply(dta,R)
        err = np.append(err, sum(sum(R)))

        print('Error: %.3e'%err[len(err)-1])

    err = err[1:]
    return u0, u, err, V, E, BE, IE