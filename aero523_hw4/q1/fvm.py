import numpy as np

def getglob(u):
    a = max(u)
    return a

def flux(ul, ur, a):
    Fhat = 1/2*ul**2
    return Fhat

def solve(x, u0, T, CFL):
    a = getglob(u0)
    dx = x[1] - x[0]
    dt = CFL*dx/a; Nt = int(np.ceil(T/dt)); dt = T/Nt
    Ne = u0.size

    u = u0.copy(); R = u.copy()
    for n in range(Nt):
        R *= 0
        for j in range(Ne+1):
            ul = u[j-1] 
            ur = u[j  ] if (j < Ne) else u[0]
            Fhat = flux(ul,ur,a)
            if (j > 0 ): R[j-1] += Fhat
            if (j < Ne): R[j  ] -= Fhat
        u -= dt/dx * R
    return u