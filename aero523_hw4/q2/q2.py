import numpy as np
import matplotlib.pyplot as plt 
from fvm import solve

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def main(save):

    State_Verify = np.transpose(np.array([1, 0.5, .25]))
    x, uRus, uHLLE, rRus, rHLLE = solve(State_Verify, State_Verify, 100, 0.2, True)
    nums  = np.arange(rRus.shape[1])

    plt.figure(figsize=(15,5))
    plt.plot(nums, rRus[0,:], linestyle='--',lw=4, label=r'HLLE - $\rho$')
    plt.plot(nums, rRus[1,:], linestyle='--',lw=4, label=r'HLLE - $\rho u$')
    plt.plot(nums, rRus[2,:], linestyle='--',lw=4, label=r'HLLE - $\rho E$')
    plt.plot(nums, rHLLE[0,:], linestyle=':', lw=4 , label=r'Rusanov - $\rho$')
    plt.plot(nums, rHLLE[1,:], linestyle=':', lw=4 , label=r'Rusanov - $\rho u$')
    plt.plot(nums, rHLLE[2,:], linestyle=':', lw=4 , label=r'Rusanov - $\rho E$')
    plt.xlabel(r'Iterations of $\Delta t$', fontsize=16)
    plt.ylabel(r'$L_2$ Residual of the norm', fontsize=16)
    plt.legend(loc='lower left', fontsize=14, ncol=6)
    if save: plt.savefig('verification.pdf', bbox_inches='tight')
    plt.show()


    ul = np.transpose(np.array([1, 0, 2.5]))
    ur = np.transpose(np.array([0.125, 0, 0.25]))
    for N in np.array([50, 100, 200]):
        x, uRus, uHLLE, rRus, rHLLE = solve(ul, ur, N, 0.2, False)

        plt.figure(figsize=(15,5))
        plt.subplot(1,3,1)
        plt.plot(x, uHLLE[0,:], lw=4, label=r'HLLE')
        plt.plot(x, uRus[0,:], linestyle='--', lw=4 , label=r'Rusanov')
        plt.xlabel(r'Distance along tube, x', fontsize=16)
        plt.ylabel(r'Density, $\rho$', fontsize=16)
        plt.legend(loc='lower left', fontsize=16)

        plt.subplot(1,3,2)
        plt.plot(x, uHLLE[1,:], lw=4, label=r'HLLE')
        plt.plot(x, uRus[1,:], linestyle='--',lw=4 , label=r'Rusanov')
        plt.xlabel(r'Distance along tube, x', fontsize=16)
        plt.ylabel(r'Momentum, $\rho u$', fontsize=16)

            
        plt.subplot(1,3,3)
        plt.plot(x, uHLLE[1,:]/uHLLE[0,:], lw=4, label=r'HLLE')
        plt.plot(x, uRus[1,:]/uRus[0,:], linestyle='--',lw=4 , label=r'Rusanov')
        plt.xlabel(r'Distance along tube, x', fontsize=16)
        plt.ylabel(r'Velocity, $ u$', fontsize=16)

        save_name = 'n' + str.format('{0:.0f}', N) + '.pdf'
        if save: plt.savefig(save_name, bbox_inches='tight')
        plt.show()
if __name__ == "__main__":
    main(True)