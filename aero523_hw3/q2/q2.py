from matplotlib import pyplot as plt
import numpy as np
import math 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



def main():    
    
    #make_plot(50, 40, 'BW_sig1', False) # num = 50, Nt = 50 -> sig = 1
    #make_plot(99, 50, 'BW_sig2', False) # Nt = 50, num = 2*Nt-1 -> sig = 2
    make_plot(150, 50, 'BW_sig2', False) # Nt = 50, num = 2*Nt-1 -> sig = 2
    #make_plot(40, 80, 'BW_stable', False)

def make_plot(num, Nt, savename, savedat):
    L = 2; a = 0.5; T = L/a
    xlin = np.linspace(0, 2, num, endpoint=True)
    us = np.zeros((Nt, num))
    for i in range(num):
        us[0,i] = math.exp(-100*(xlin[i]/L - 0.5)**2)  
    
    dx = xlin[1]; dt = T/(Nt-1)
    sig = a*dt/dx
    print('Sig = ', sig)
   

    plt.figure(figsize=(9,5))
    plt.plot(np.NaN, np.NaN, '-', color='none', label = r'$\sigma$ = ' + str.format('{0:.2f}', sig))
    for n in range(Nt):
        if (n+1)%10 == 0 or n == 0:
            plot_label = r'$t$ = ' + str.format('{0:.2f}', n*dt)
            plt.plot(xlin, us[n,:], lw = 2, label = plot_label)
        if n < Nt-1:
            for j in range(2, num):
                uj = us[n,j]
                ujm1 = us[n,j-1]
                ujm2 = us[n,j-2]
                us[n+1,j] = us[n,j] - sig/2*(3*uj - ujm1 + ujm2) + sig**2/2*(ujm2 - 2*ujm1 + uj)
    plt.xlabel(r'Location along $x$', fontsize = 16)
    plt.ylabel(r'BW Approximation, $u$', fontsize = 16)
    plt.legend(loc = 'center left', fontsize = 18)
    if savedat: plt.savefig(savename + '.pdf', bbox_inches = 'tight')
    plt.show()

if __name__ == "__main__":
    main()