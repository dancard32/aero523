def directsol(p):
    import numpy as np
    from scipy import sparse
    from scipy.sparse import linalg

    N = 2**(p + 3)  # P-value scaling
    h = float(2/N)  # Step size
    xlin = np.linspace(-1, 1, N+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values
    
    A = sparse.lil_matrix(((N+1)**2,(N+1)**2))      # Pre-allocate sparse Matrix
    q = np.zeros((N+1)**2)                          # Pre-allocate q vector
    for iy in range(N+1):
        for ix in range(N+1):
            i = iy*(N+1) + ix                       # Iteration index
            iL = i - 1; iR = i + 1                  # Left/Right indices
            iD = i - (N+1); iU = i + (N+1)          # Top/Bottom indices

            if ylin[iy,ix] == 1:     # Top Boundary
                q[i] = 1
                A[i,i] = 1
            elif ylin[iy,ix] == -1:  # Bottom Boundary
                q[i] = -1
                A[i,i] = 1
            elif abs(xlin[iy,ix]) == 1:  # Left/Right Boundary
                q[i] = ylin[iy,ix]
                A[i,i] = 1
            elif abs(xlin[iy,ix]) <= 0.25 and abs(ylin[iy,ix]) <= 0.25: # Interior Boundary
                q[i] = 0
                A[i,i] = 1
            else:   # Interior Domain
                q[i] = 0
                A[i,i] = -4*h**-2
                A[i,iL] = h**-2
                A[i,iR] = h**-2
                A[i,iD] = h**-2
                A[i,iU] = h**-2
    phiv = linalg.spsolve(sparse.csr_matrix(A), q)  # Solve for Phi
    phi = np.reshape(phiv, (N+1, N+1))              # Re-shape for plotting

    return xlin, ylin, phi
