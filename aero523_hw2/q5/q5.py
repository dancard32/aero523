import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def  advec_diff(N, nu):
    L = 1
    x = np.linspace(0, L, N+1)
    dx = float(L)/N
    data = np.ones((3, N-1))
    data[0,:] *= -(dx + 2*nu); data[1,:] *= 4*nu; data[2,:] *= (dx - 2*nu)
    diags = np.array([-1, 0, 1])
    
    A = sparse.spdiags(data, diags, N-1, N-1, 'csr')
    b = np.zeros(N+1)

    q = 2*dx**2*np.sin(np.pi * x)
    b = q[1:N]
    b[0] += 0
    b[N-2] += 0

    Tt = linalg.spsolve(A, b)
    T = np.zeros(N+1)
    T[0] = 0; T[1:N] = Tt; T[N] = 0
    return x, T

x, T = advec_diff(16, 0.1)
x2, T2 = advec_diff(16, 1)

plt.figure(figsize=(8,4))
plt.plot(x, T, color = 'black', lw = 2, label=r'$\nu = 0.1$')
plt.plot(x2, T2, color = 'gray', lw = 2, label=r'$\nu = 1$')
plt.xlabel(r'X-Axis', fontsize = 16)
plt.ylabel(r'Approximated $u$', fontsize = 16)
plt.legend(fontsize = 18)
plt.savefig('nu_values.pdf', bbox_inches = 'tight')
plt.show()



Ns = [4, 8, 16, 32]
plt.figure(figsize=(8,4))
for i in range(len(Ns)):
    x, T = advec_diff(Ns[i], 0.1)
    plot_label = r'$N$ = ' + str(Ns[i])
    plt.plot(x, T, lw = 2, label = plot_label)
plt.xlabel(r'X-Axis', fontsize = 16)
plt.ylabel(r'Approximated $u$', fontsize = 16)
plt.legend(fontsize = 18)
plt.savefig('varying_Ns.pdf', bbox_inches = 'tight')
plt.show()