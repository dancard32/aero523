from matplotlib import pyplot as plt
import numpy as np
import math 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def main():    
    
    nt = 1000    
    xints = np.arange(100, 2010, step=10, dtype=int)
    err_dat = np.zeros(xints.shape[0])
    k = 0
    for nx in xints:
        print(nx)
        us, xlin = BW_method(nx, nt)    # Call BW function

        u0 = us[0,:]                # First value from u
        ul = us[us.shape[0]-1,:]    # Last value from u
        diff = abs(ul - u0)         # Take difference
        err_dat[k] = math.sqrt(np.dot(diff, diff)/(nx+1)) # L2 norm
        
        k += 1
    
    endval = (len(err_dat) - 1)//2 #Take the approximate midpoint
    rate = math.log10(err_dat[endval]/err_dat[0])/math.log10(xints[endval]/xints[0])
    print('Nt-', nt, ': converge rate = ', abs(rate))
    
    plt.figure(figsize=(11,6))
    plt.plot(xints, err_dat, lw=2, label = r'$N_t$ = ' + str.format('{0:.0f}', nt))
    plt.xlabel(r'$N_x$ Intervals', fontsize = 20)
    plt.ylabel(r'$L_2$ Residual Error Norm', fontsize = 20)
    plt.yscale('log')
    plt.xscale('log')
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.ylim(10**-6.5, 10**(0.5))
    plt.legend(loc = 'upper center', fontsize = 18, ncol=4)
    #plt.savefig('BW_convergence.pdf', bbox_inches = 'tight')
    plt.show()
    
    
    """
    plt.figure(figsize=(11,6))
    us, xlin = BW_method(100, 1000)
    plt.plot(xlin, us[0, :], lw=2, label='u(x,0)')
    plt.plot(xlin, us[len(us)-1, :], lw=2, label='u(x,T)')
    plt.legend(loc = 'upper right')
    plt.show()
    print(us[0,:])
    print(us[us.shape[0]-1,:])
    """
def BW_method(nx, nt):
    L = 2; a = 0.5; T = L/a
    xlin = np.linspace(0, 2, nx+1, endpoint=True)
    us = np.zeros((nt+1, nx+1))
    
    us[0,:] = np.exp(-100*(xlin/L - 0.5)**2)   
    dx = xlin[1]; dt = T/(nt-1)
    sig = a*dt/dx
    
    for n in range(0, nt):
        for j in range(0, nx+1):
            ujm1 = us[n,(j-1)%(nx+1)]
            ujm2 = us[n,(j-2)%(nx+1)]
            uj = us[n,j]

            us[n+1,j] = us[n,j] - sig/2*(3*uj - 4*ujm1 + ujm2) + sig**2/2*(ujm2 - 2*ujm1 + uj)

    return us, xlin

if __name__ == "__main__":
    main()