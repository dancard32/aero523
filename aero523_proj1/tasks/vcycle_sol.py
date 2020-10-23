import numpy as np
import math 
from direct_sol import directsol

def residual(U, F):
    N = U.shape[0] - 1; h = 2.0/N
    R = U.copy()

    for iy in range(1, N):
        for ix in range(1, N):
            R[iy, ix] = F[iy, ix] - (U[iy+1, ix] + U[iy-1, ix] + U[iy, ix+1] + U[iy, ix-1] - 4*U[iy, ix])*h**-2
    return R
def restrict(r, N):
    rc = np.zeros([int(N/2)+1, int(N/2)+1])
    
    ix = 0; iy = 0
    for j in range(1, N+1, 2):
        for i in range(1, N+1, 2):            
            rc[iy,ix] =  1/8.*(r[j+1,i] + r[j-1,i] + r[j,i+1] + r[j,i-1]) + 1/16.*(r[j+1,i+1] + r[j-1,i-1] + r[j-1,i+1] + r[j+1,i-1]) + 1/4.*r[j,i]
            ix += 1
            if ix >= int(N/2):
                ix = 0
        iy += 1
    return rc
def prolongate(e2h, N):
    eh = np.zeros([int(2*N)+1, int(2*N)+1])
    xlin = np.linspace(-1, 1, int(2*N)+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values

    for j in range(1, int(2*N), 2):
        for i in range(1, int(2*N), 2):
                if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                    pass
                else:   # Interior Domain
                    # Self Weight
                    eh[j, i] =  e2h[math.floor(j/2.), math.floor(i/2.)]

                    # Up/Down Nodes
                    eh[j+1, i] += 0.5*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j-1, i] += 0.5*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j, i+1] += 0.5*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j, i-1] += 0.5*e2h[math.floor(j/2.), math.floor(i/2.)]

                    # Corner Nodes
                    eh[j+1, i+1] += 0.25*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j-1, i-1] += 0.25*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j-1, i+1] += 0.25*e2h[math.floor(j/2.), math.floor(i/2.)]
                    eh[j+1, i-1] += 0.25*e2h[math.floor(j/2.), math.floor(i/2.)]
    return eh

def multigrid(U, F, p, pmax, viter, nu1, nu2, nuc):
    l2err = np.zeros(viter)
    fh = F

    print("V-Cycle Method(p=", pmax, ")\n---------------------")
    for k in range(viter):
        print("Iteration: ", k)
        
        N = U.shape[0] - 1
        griditer = 0
        fmat = np.zeros([U.shape[0], U.shape[0], pmax-p+1])

        # Sweep Down
        utemp = smooth(U, fh, nu1)
        while N > 2**(p + 3):
            rh = residual(utemp, fh)
            f2h = restrict(rh, N)
            
            griditer += 1
            N = int(N/2)
            fmat[0:(N+1), 0:(N+1), griditer] = f2h
            utemp = smooth(np.zeros([N+1, N+1]), f2h, nu1)

        # Coarsest Mesh
        utemp = smooth(np.zeros([N+1, N+1]), fmat[0:(N+1), 0:(N+1), griditer], nuc)

        # Sweep Up
        while N <= 2**(pmax + 2):
            utemp = prolongate(utemp, N)

            N = int(2*N)
            griditer -= 1
            utemp = smooth(utemp, fmat[0:(N+1), 0:(N+1), griditer], nu2)

        U += utemp
        resid = residual(U, F)
        
        for j in range(N + 1):
            for i in range(N + 1):
                l2err[k] += resid[j,i]**2
        l2err[k] = np.sqrt(l2err[k]/(N + 1)**2)

    return U, l2err
def vcyclesol(p, pmax, viter, nu1, nu2, nuc):
    N = 2**(pmax + 3)
    U = np.zeros([N+1, N+1])
    F = U.copy()
    U[:, 0] = np.flip(np.linspace(-1, 1, N+1))
    U[:, N] = np.flip(np.linspace(-1, 1, N+1))
    U[0,:] = 1; U[N,:] = -1    

    U, l2err = multigrid(U, F, p, pmax, viter, nu1, nu2, nuc)

    return U, l2err
def smooth(U, F, nu):
    N = U.shape[0]-1; h = 2.0/N    
    omega = 1.5
        
    xlin = np.linspace(-1, 1, N+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values
    for k in range(nu):
        for iy in range(1, N):   # Red Nodes --------------------------------------------
            if iy%2 == 0:
                for ix in range(1, N, 2):
                    if abs(xlin[iy,ix]) <= 0.25 and abs(ylin[iy,ix]) <= 0.25: # Interior Boundary
                        U[iy, ix] = 0
                    else:   # Interior Domain
                        unew = 0.25*(U[iy+1, ix] + U[iy-1, ix] + U[iy, ix-1] + U[iy, ix+1] - F[iy, ix]*h**2)
                        U[iy, ix] = U[iy, ix]*(1.0 - omega) + unew*omega
            else:
                for ix in range(2, N-1, 2):
                    if abs(xlin[iy,ix]) <= 0.25 and abs(ylin[iy,ix]) <= 0.25: # Interior Boundary
                        U[iy, ix] = 0
                    else:   # Interior Domain
                        unew = 0.25*(U[iy+1, ix] + U[iy-1, ix] + U[iy, ix-1] + U[iy, ix+1] - F[iy, ix]*h**2)
                        U[iy, ix] = U[iy, ix]*(1.0 - omega) + unew*omega
        for iy in range(1, N):   # Black Nodes --------------------------------------------
            if iy%2 == 0:
                for ix in range(2, N-1, 2):
                    if abs(xlin[iy,ix]) <= 0.25 and abs(ylin[iy,ix]) <= 0.25: # Interior Boundary
                        U[iy, ix] = 0
                    else:   # Interior Domain
                        unew = 0.25*(U[iy+1, ix] + U[iy-1, ix] + U[iy, ix-1] + U[iy, ix+1] - F[iy, ix]*h**2)
                        U[iy, ix] = U[iy, ix]*(1.0 - omega) + unew*omega
            else:
                for ix in range(1, N, 2):
                    if abs(xlin[iy,ix]) <= 0.25 and abs(ylin[iy,ix]) <= 0.25: # Interior Boundary
                        U[iy, ix] = 0
                    else:   # Interior Domain
                        unew = 0.25*(U[iy+1, ix] + U[iy-1, ix] + U[iy, ix-1] + U[iy, ix+1] - F[iy, ix]*h**2)
                        U[iy, ix] = U[iy, ix]*(1.0 - omega) + unew*omega
    return U