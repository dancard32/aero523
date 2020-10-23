def gausssol(num_iters, p, omega):
    import numpy as np
    import math
    # Variables for Jacobi Iteration
    N = 2**(p + 3)  # P-value scaling
    h = float(2/N)  # Step size
    xlin = np.linspace(-1, 1, N+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values
    resid_mat = np.zeros((N+1, N+1))
    resid = np.zeros(num_iters)

    # Boundary Conditions / Intial Guess
    U = np.zeros((N+1, N+1))
    U[0,:] = 1; U[N,:] = -1; U[:,0] = ylin[:,0]; U[:,N] = ylin[:,0]
    for k in range(num_iters):
        for j in range(1, N):   # Red Nodes --------------------------------------------
            if j%2 == 0:
                for i in range(1, N, 2):
                    if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                        resid_mat[j,i] = 0
                    else:   # Interior Domain
                        u_unew = (U[j+1, i] + U[j-1, i] + U[j,i+1] + U[j,i-1])/4
                        rij = 0 - h**-2*(4*u_unew - 4*U[j,i])
                        resid_mat[j,i] = rij
                        U[j,i] = U[j,i]*(1.-omega) + u_unew*omega
            else:
                for i in range(2, N-1, 2):
                    if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                        resid_mat[j,i] = 0
                    else:   # Interior Domain
                        u_unew = (U[j+1, i] + U[j-1, i] + U[j,i+1] + U[j,i-1])/4
                        rij = 0 - h**-2*(4*u_unew - 4*U[j,i])
                        resid_mat[j,i] = rij
                        U[j,i] = U[j,i]*(1.-omega) + u_unew*omega
        for j in range(1, N):   # Black Nodes --------------------------------------------
            if j%2 == 0:
                for i in range(2, N-1, 2):
                    if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                        resid_mat[j,i] = 0
                    else:   # Interior Domain
                        u_unew = (U[j+1, i] + U[j-1, i] + U[j,i+1] + U[j,i-1])/4
                        rij = 0 - h**-2*(4*u_unew - 4*U[j,i])
                        resid_mat[j,i] = rij
                        U[j,i] = U[j,i]*(1.-omega) + u_unew*omega
            else:
                for i in range(1, N, 2):
                    if abs(xlin[j,i]) <= 0.25 and abs(ylin[j,i]) <= 0.25: # Interior Boundary
                        resid_mat[j,i] = 0
                    else:   # Interior Domain
                        u_unew = (U[j+1, i] + U[j-1, i] + U[j,i+1] + U[j,i-1])/4
                        rij = 0 - h**-2*(4*u_unew - 4*U[j,i])
                        resid_mat[j,i] = rij
                        U[j,i] = U[j,i]*(1.-omega) + u_unew*omega
        for i in range(N+1):
            for j in range(N+1):
                resid[k] += resid_mat[j,i]**2
        resid[k] = math.sqrt(resid[k]/(N+1)**2)
    return resid