from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import math 
from direct_sol import directsol
from calc_coefs import calc_cp, calc_cl
from jacobi_sol import jacobisol
from gauss_sol import gausssol
from vcycle_sol import vcyclesol

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def export_phi(phi):
    f = open('9by9_mat',"w")  # Filename
    output = ''
    for j in range(9):
        for i in range(9):
            output += str.format('{0:.4f}',phi[j,i]) + r'& '
        if i == 8:
            output +=  r'\\' # Output results to LaTeX environment
    f.write(output)
    f.close()

def export_cl(cl):
    f = open('cl_vals',"w")  # Filename
    output = ''
    for i in range(6):
        output += r'$p = $ ' + str.format('{0:.0f}',i)+ r'& ' + str.format('{0:.7f}',cl[i]) + r'\\'
    f.write(output)
    f.close()

def gen_grids(p):
    N = 2**(p + 3)  # P-value scaling
    xlin = np.linspace(-1, 1, N+1, endpoint=True)   # linspace over domain
    xlin, ylin = np.meshgrid(xlin, np.flip(xlin))   # Meshgrid values

    return xlin, ylin

def run_q1():
    xlin, ylin, phi = directsol(0)
    export_phi(phi)

def run_q2():
    plt.figure(figsize=(8,4))
    for i in np.array([0, 2, 4]):
        xlin, ylin, phi = directsol(i)
        cp = calc_cp(ylin, phi)

        plot_label = r'$p$ = ' + str(i)
        plt.plot(xlin[0,:], -cp, lw = 2, label = plot_label)
    plt.xlabel(r'Location along bottom wall', fontsize = 16)
    plt.ylabel(r'$-c_p(x)$', fontsize = 16)
    plt.legend(fontsize = 18)
    plt.savefig('figs/cp_runs.pdf', bbox_inches = 'tight')
    plt.show()

    xlin, ylin, phi = directsol(7)
    cl_exact = calc_cl(ylin, calc_cp(ylin, phi))
    f = open('cl_exact',"w") 
    output = str.format('{0:.7f}',cl_exact)
    f.write(output)
    f.close()

    cl_vals = np.zeros(6); h_vals = np.zeros(6); cl_err = np.zeros(6)
    for i in range(6):
        xlin, ylin, phi = directsol(i)
        cl = calc_cl(ylin, calc_cp(ylin, phi))
        cl_vals[i] = cl
        h_vals[i] = ylin[0,0] - ylin[1,0]
        cl_err[i] = abs(cl - cl_exact)
    export_cl(cl_vals)
    
    rate = math.log10(cl_err[4]/cl_err[5])/math.log10(h_vals[4]/h_vals[5])
    plot_label = r'Convergence = $\mathcal{O}($' + str.format('{0:.4f}',rate) + R')'
    plt.figure(figsize=(8,4))
    plt.plot(range(6), cl_err, color = 'black', marker = 'o', lw = 2, label = plot_label)
    plt.xlabel(r'$p$ values', fontsize = 16)
    plt.ylabel(r'Error, $|| c_l(p) - c_{l, exact}||$', fontsize = 16)
    plt.yscale('log')
    plt.legend(fontsize = 18)
    plt.savefig('figs/cl_err.pdf', bbox_inches = 'tight')
    plt.show()

def run_q3():
    plt.figure(figsize=(8,4))
    plt.plot(range(2500), jacobisol(0.3), color = 'black', lw = 2, label = r'$\omega = 0.3$'); print('Omega = 0.3 - Done')
    plt.plot(range(2500), jacobisol(0.6), color = 'blue', lw = 2, label = r'$\omega = 0.6$'); print('Omega = 0.6 - Done')
    plt.plot(range(2500), jacobisol(1.0), color = 'gray', lw = 2, label = r'$\omega = 1.0$'); print('Omega = 1.0 - Done')
    plt.xlabel(r'Iteration Number', fontsize = 16)
    plt.ylabel(r'$L_2$ Residual Norm Error', fontsize = 16)
    plt.yscale('log')
    plt.legend(fontsize = 18)
    plt.savefig('figs/jacobi_l2.pdf', bbox_inches = 'tight')
    plt.show()

    plt.figure(figsize=(8,4))
    plt.plot(range(2500), gausssol(2500, 3, 0.5), color = 'black', lw = 2, label = r'$\omega = 0.5$'); print('Omega = 0.5 - Done')
    plt.plot(range(2500), gausssol(2500, 3, 1.0), color = 'blue', lw = 2, label = r'$\omega = 1.0$'); print('Omega = 1.0 - Done')
    plt.plot(range(2500), gausssol(2500, 3, 1.5), color = 'gray', lw = 2, label = r'$\omega = 1.5$'); print('Omega = 1.5 - Done')
    plt.xlabel(r'Iteration Number', fontsize = 16)
    plt.ylabel(r'$L_2$ Residual Norm Error', fontsize = 16)
    plt.yscale('log')
    plt.legend(fontsize = 18)
    plt.savefig('figs/gauss_l2.pdf', bbox_inches = 'tight')
    plt.show()    

def run_q4():    
    p = 1; viter = 25
    ax = plt.figure(figsize=(10,5)).gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for pval in range(2, 6):
        U, l2_2 = vcyclesol(p, pval, viter, 2, 2, 50)
        plotlabel = r'$p$ = ' + str(pval)
        plt.plot(range(1, viter+1), l2_2, lw = 2, label = plotlabel)
    plt.xlabel(r'Multi-grid Cycles', fontsize = 16)
    plt.ylabel(r'$L_2$ Residual Norm', fontsize = 16)
    plt.yscale('log')
    plt.legend(fontsize = 18, ncol = 4)
    plt.savefig('figs/vcyc_l2_err.pdf', bbox_inches = 'tight')
    plt.show()
    
    
    p = 0; pmax = 4; viter = 50
    nu1 = 10; nu2 = 10; nuc = 1000
    U, resid_norm = vcyclesol(p, pmax, viter, nu1, nu2, nuc)
    workunits = np.zeros(viter)
    for l in range(viter-1):
        workunits[l+1] = workunits[l] + 2**(pmax*l/viter + 3)/(nu1 + nu2)
        
    plt.figure(figsize=(8,4))
    plt.plot(workunits, resid_norm, color = 'black', lw = 2, label = r'Multigrid')
    plt.plot(range(100), gausssol(100, 3, 1.5), color = 'gray', lw = 2, label = r'Gauss-Seidel')
    plt.xlabel(r'Work Units', fontsize = 16)
    plt.ylabel(r'$L_2$ Residual Norm', fontsize = 16)
    plt.yscale('log')
    plt.legend(fontsize = 18)
    plt.savefig('figs/vcyc_gs.pdf', bbox_inches = 'tight')
    plt.show()   
def main():    
    run_q1()
    run_q2()
    run_q3()
    run_q4()
if __name__ == "__main__":
    main()