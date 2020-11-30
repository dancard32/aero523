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
    
def calcATPR(u0, u, alpha, V, BE):
    gam = 1.4

    Pinf = (gam-1)*(u0[0,3]-0.5*u0[0,0]*((u0[1,0]/u0[0,0])**2 + (u0[2,0]/u0[0,0])**2))
    Ptinf = Pinf*(1 + 0.5*(gam-1)*(2.2)**2)*(gam/(gam-1))

    ATPR = 0; d = 0
    for i in range(BE.shape[0]):
        n1, n2, e1, bgroup = BE[i,:]
        xl = V[n1,:]; xr = V[n2,:]
        uedge = u[e1,:]

        dy = xr[1] - xl[1]
        
        if bgroup == 1: # Exit
            uvel = uedge[1]/uedge[0]; vvel = uedge[2]/uedge[0]
            q = np.sqrt(uvel**2 + vvel**2)
            P = (gam-1)*(uedge[3]-0.5*uedge[0]*q**2)  
            c = np.sqrt(gam*P/uedge[0])
            mach = q/c
            Pt = P*(1 + 0.5*(gam-1)*mach**2)**(gam/(gam-1))

            d += dy
            ATPR += Pt*dy/Ptinf
            
    ATPR *= 1/d
    return ATPR
                
def solve(alpha, mesh):
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    u0 = getIC(alpha, E.shape[0]); u = u0.copy(); ATPR = np.array([calcATPR(u0,u,1,V,BE)])
    R = np.zeros((E.shape[0], 4)); dta = R.copy(); err = np.array([1]); itr = 0

    while err[err.shape[0]-1] > 10**(-5):
    #for k in range(50):
        R *= 0; dta *= 0
        for i in range(IE.shape[0]):
            n1, n2, e1, e2 = IE[i,:]
            xl = V[n1,:]; xr = V[n2,:]
            ul = u[e1,:]; ur = u[e2,:]

            dx = xr - xl; deltal = LA.norm(dx)
            nhat = np.array([dx[1], -dx[0]])/deltal
            F, FL, FR, ls = RoeFlux(ul, ur, nhat)
            R[e1,:] += F*deltal; R[e2,:] -= F*deltal
            dta[e1,:] += ls*deltal; dta[e2,:] += ls*deltal

        for i in range(BE.shape[0]):
            n1, n2, e1, bgroup = BE[i,:]
            xl = V[n1,:]; xr = V[n2,:]
            uedge = u[e1,:]

            dx = xr - xl; deltal = LA.norm(dx)
            nhat = np.array([dx[1], -dx[0]])/deltal

            if bgroup == 0: # Engine - Invscid
                vp = np.array([uedge[1], uedge[2]])/uedge[0]
                vb = vp - np.dot(vp, nhat)*nhat
                pb = 0.4*(uedge[3] - 0.5*uedge[0]*(vb[0]**2 + vb[1]**2))
                ignore, FL, FR, ls = RoeFlux(uedge, u0[0,:], nhat)

                F = pb*np.array([0, nhat[0], nhat[1], 0])
            elif bgroup == 1 or bgroup == 2: # Exit/Outflow - Supersonic Outflow
                F, FL, FR, ls = RoeFlux(uedge, uedge, nhat)
            elif bgroup == 3: # Inflow
                F, FL, FR, ls = RoeFlux(uedge, u0[0,:], nhat)
            
            R[e1,:] += F*deltal
            dta[e1,:] += ls*deltal
        
        dta = 2/dta
        u -= np.multiply(dta, R)
        err = np.append(err, sum(sum(abs(R))))

        ATPR = np.append(ATPR, calcATPR(u0,u,1,V,BE))
        print('Iteration: %3d,\tError: %.3e, ATPR: %.3f'%(itr, err[err.shape[0]-1], ATPR[ATPR.shape[0]-2])); itr += 1

    return u, err[1:], ATPR, V, E, BE, IE