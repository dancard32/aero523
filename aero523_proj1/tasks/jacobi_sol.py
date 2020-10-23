def jacobisol(omega):
    import numpy as np
    import math
    
    p = 3
    N = 2**(p + 3)  # P-value scaling
    h = float(2/N)  # Step size
    xlin = np.linspace(-1, 1, N+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values

    # Variables for Jacobi Iteration
    resid_mat = np.zeros((N+1, N+1))
    num_iters = 2500
    resid = np.zeros(num_iters)

    # Boundary Conditions / Intial Guess
    U = np.zeros((N+1, N+1))
    U[0,:] = 1; U[N,:] = -1
    U[:,0] = ylin[:,0]; U[:,N] = ylin[:,0]
    u_temp = U.copy()

    for k in range(num_iters):
        for j in range(1, N):
            for i in range(1, N):
                if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                    u_temp[j,i] = 0
                    resid_mat[j,i] = 0
                else:   # Interior Domain
                    u_temp[j,i] = (U[j+1, i] + U[j-1, i] + U[j,i+1] + U[j,i-1])/4
                    rij = 0 - h**-2*(4*u_temp[j,i] - 4*U[j,i])
                    resid_mat[j,i] = rij

        U[1:N, 1:N] = U[1:N, 1:N]*(1.-omega) + u_temp[1:N, 1:N]*omega
        for i in range(N+1):
            for j in range(N+1):
                resid[k] += resid_mat[j,i]**2
        resid[k] = math.sqrt(resid[k]/(N+1)**2)

    return resid
