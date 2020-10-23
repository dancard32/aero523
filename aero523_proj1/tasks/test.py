import numpy as np
from scipy import sparse
from scipy.sparse import linalg
from matplotlib import pyplot as plt
import array


test = np.arange(16).reshape(4,4)
#test = np.delete(test, 2, 1)
#test = np.delete(test, 2, 0)

print(test)
test = np.hstack((test[:, 0:2], test[:, 3:4]))


print('Split')
print(test)
test = np.vstack((test[0:2,:], test[3:4, :]))

print('Split')
print(test)

    p = 1; pmax = 2; viter = 1; nu1 = 10; nu2 = 10; nuc = 50
    U, l2 = vcyclesol(p, pmax, viter, nu1, nu2, nuc)

    
    xlin, ylin = gen_grids(pmax)

    gradients = np.linspace(-1, 1, 50)
    plt.figure()
    plt.contourf(xlin, ylin, U, gradients, cmap = 'jet')
    plt.fill_between(np.array([-0.25, 0.25]), np.array([-0.25, -0.25]), np.array([0.25, 0.25]), color = 'black')
    plt.axis('equal')
    plt.colorbar()
    plt.show()