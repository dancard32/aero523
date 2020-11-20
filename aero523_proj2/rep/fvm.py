import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from flux import RoeFlux
from readgri import readgri, writegri

def getIC(alpha, Ne):
    alpha = np.deg2rad(alpha); Minf = 2.2; gam = 1.4
    uinf = np.array([1, Minf*np.cos(alpha), Minf*np.sin(alpha), 1/(gam*(gam-1)) + Minf**2/2])

    u0 = np.zeros((Ne, 4))
    for i in range(4):
        u0[:,i] = uinf[i]
    u0[abs(u0) < 10**-10]

    return u0

def solve():
    mesh = readgri('mesh0.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    # Get the initial state
    u0 = getIC(1, E.shape[0]); u = u0.copy()

    # Initialize the residual matrix, dti/Ai, and error
    R = np.zeros((E.shape[0], 4)); err = np.array([]); dta = np.zeros(E.shape[0])
    
    plt.figure()
    print('Residual Norm\n' + 15*'-'); 
    # Iterate while out of tolerance
    for itr in range(1000):
        R *= 0; dta *= 0    # Re-initialize

        # Loop over Interior Edges
        for i in range(IE.shape[0]):
            n1, n2, e1, e2 = IE[i,:]
            x1 = V[n1,:]; x2 = V[n2,:]
            u1 = u[e1,:]; u2 = u[e2,:]

            dx = x2 - x1; deltal = np.sqrt(dx[0]**2 + dx[1]**2)
            nhat = np.array([-dx[1],dx[0]])/deltal
            
            # Determine the flux
            F, ls = RoeFlux(u1, u2, nhat, False)
            R[e1,:] += F*deltal; R[e2,:] -= F*deltal
            dta[e1] += ls*deltal; dta[e2] += ls*deltal

        # Implement Boundaries - FREE STREAM TEST
        for i in range(BE.shape[0]):
            n1, n2, elem, bgroup = BE[i,:]
            x1 = V[n1,:]; x2 = V[n2,:]
            uedge = u[elem,:]

            dx = x2 - x1; deltal = np.sqrt(dx[0]**2 + dx[1]**2)
            nhat = np.array([dx[1],-dx[0]])/deltal

            # Engine Boundary Condition
            if bgroup == 0:
                
                # Inviscid Boundary Conditions
                vplus = np.array([uedge[1], uedge[2]])/uedge[0]
                vb = vplus - np.dot(vplus, nhat)*nhat
                pb = (1.4-1)*(uedge[3] - 0.5*uedge[0]*(vb[0]**2 + vb[1]**2))

                # Call roe flux for eigenvalues of wave propagation (ls)
                ignore, ls = RoeFlux(uedge, uedge, nhat, False) 

                # Determine the flux
                F = np.array([0, pb*nhat[0], pb*nhat[1], 0])

                F, ls = RoeFlux(uedge,  u0[0,:], nhat, False)  # <- Un-comment for free-stream test
                R[elem,:] -= F*deltal

            # Exit and Outflow Boundary Condition
            elif bgroup == 1 or bgroup == 2:

                # Down-stream conditions are equal
                F, ls = RoeFlux(uedge, uedge, nhat, False)
                R[elem,:] -= F*deltal
            
            # Inflow Boundary Condition
            elif bgroup == 3:
                
                # Set to the initial state for inflow
                F, ls = RoeFlux(uedge, u0[0,:], nhat, False)
                R[elem,:] -= F*deltal
            
            dta[elem] += ls*deltal

        dta = 2*0.5/dta
        for i in range(E.shape[0]):
            u[i] -= dta[i]*R[i,:]
        err = np.append(err, sum(sum(abs(R))))


        P = (1.4 - 1)*(u[:,3] - 0.5*u[:,0]*((u[:,1]/u[:,0])**2 + (u[:,2]/u[:,0])**2))
            
        plt.cla()
        plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=P, cmap='jet', shading='flat')
        plt.axis('equal'); plt.axis('off')
        plt.draw()
        plt.pause(0.1)


        print('Iteration %3d, Error: %.3e'%(itr, err[len(err)-1]))

    return u0, u, err, V, E, BE, IE
