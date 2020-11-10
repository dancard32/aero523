import matplotlib.pyplot as plt
import numpy as np
import math
from fvm import solve, getglob, flux

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def q1a():
    
    plt.figure(figsize=(6,6))
    for t in np.array([0.5, 1.0, 1.5]):
        x0, u0 = analytical(1000, t)
        plot_label = r't=' + str(t)
        plt.plot(x0, u0, lw = 2, label=plot_label)
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'Analytical state, $u(x,t)$', fontsize=16)
    plt.legend(loc='upper left', fontsize=16)
    plt.savefig('analytical.pdf', bbox_inches='tight')
    plt.show()
    

    
    plt.figure(figsize=(6,6))
    plt.plot(np.NaN, np.NaN, lw=2, color='g')
    plt.plot(np.NaN, np.NaN, lw=2, color='r')
    plt.plot(np.NaN, np.NaN, lw=2, color='m')
    plt.plot(np.NaN, np.NaN, lw=2, color='b')
    plt.plot(np.NaN, np.NaN, lw=2, color='k')
    characteristics('1')
    characteristics('2')
    characteristics('3')
    characteristics('4')
    t = np.linspace(1, 1.5, num=50, endpoint=True); xs = 1 + np.sqrt(2 + 2*t)
    plt.plot(xs, t, lw=2, color='k')
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'time axis', fontsize=16)
    plt.legend([r'Region 1', r'Region 2', r'Region 3', r'Region 4', r'Shock Path'], loc='center left', fontsize=16)
    plt.savefig('characteristics.pdf', bbox_inches='tight')
    plt.show()

def q1b():
    
    xex, uex05 = analytical(1000, 0.5)
    xex, uex10 = analytical(1000, 1.0)
    xex, uex15 = analytical(1000, 1.5)

    

    Nx = 128; 
    x128, u0 = getIC(Nx)
    uNx12805 = solve(x128, u0, 0.5, 0.8)
    uNx12810 = solve(x128, u0, 1.0, 0.8)
    uNx12815 = solve(x128, u0, 1.5, 0.8)

    Nx = 512; 
    x512, u0 = getIC(Nx)
    uNx51205 = solve(x512, u0, 0.5, 0.8)
    uNx51210 = solve(x512, u0, 1.0, 0.8)
    uNx51215 = solve(x512, u0, 1.5, 0.8)

    plt.figure(figsize=(8,3))
    plt.plot(xex, uex05, lw=2, color='r', label=r'$u_{exact}(x,t=0.5)$')
    plt.plot(x128, uNx12805, lw=2, color='b', label=r'$N_x = 128,\ u(x,t=0.5)$')
    plt.plot(x512, uNx51205, lw=2, color='g', label=r'$N_x = 512,\ u(x,t=0.5)$')
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'$u(x,t)$ approximation', fontsize=16)
    plt.legend(loc='upper left', fontsize=16)
    plt.savefig('t05.pdf', bbox_inches = 'tight')
    plt.show()

    plt.figure(figsize=(8,3))
    plt.plot(xex, uex10, lw=2, color='r', label=r'$u_{exact}(x,t=1.0)$')
    plt.plot(x128, uNx12810, lw=2, color='b', label=r'$N_x = 128,\ u(x,t=1.0)$')
    plt.plot(x512, uNx51210, lw=2, color='g', label=r'$N_x = 512,\ u(x,t=1.0)$')
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'$u(x,t)$ approximation', fontsize=16)
    plt.legend(loc='upper left', fontsize=16)
    plt.savefig('t10.pdf', bbox_inches = 'tight')
    plt.show()

    plt.figure(figsize=(8,3))
    plt.plot(xex, uex15, lw=2, color='r', label=r'$u_{exact}(x,t=1.5)$')
    plt.plot(x128, uNx12815, lw=2, color='b', label=r'$N_x = 128,\ u(x,t=1.5)$')
    plt.plot(x512, uNx51215, lw=2, color='g', label=r'$N_x = 512,\ u(x,t=1.5)$')
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'$u(x,t)$ approximation', fontsize=16)
    plt.legend(loc='upper left', fontsize=16)
    plt.savefig('t15.pdf', bbox_inches = 'tight')
    plt.show()

def q1c():
    nxs = np.array([128, 256, 512, 1024])
    errs = np.zeros(4); k = 0

    for nx in nxs:
        x, u0 = getIC(nx)
        u = solve(x, u0, 0.5, 0.8)
        xex, uex = analytical(nx, 0.5)

        errs[k] = l2err(u, uex); k += 1

    f = open('convergences', 'w'); output = ''
    for i in range(3):
        rate = abs(math.log10(errs[i+1]/errs[i])/math.log10(nxs[i+1]/nxs[i]))
        output += r'$N_x = $ ' + str.format('{0:.0f}',nxs[i])+ r'& ' + str.format('{0:.3f}', rate) + r'\\'
    f.write(output)
    f.close()

    rate = math.log10(errs[1]/errs[0])/math.log10(nxs[1]/nxs[0])
    plotlabel = r'rate = $\mathcal{O}$(' + str.format('{0:.3f}', abs(rate)) + ')'

    plt.figure(figsize=(8,3))
    plt.plot(nxs, errs, lw=2, color='k', label=plotlabel)
    plt.xlabel(r'$N_x$ intervals', fontsize=16)
    plt.ylabel(r'$L_2$ solution error at $t=0.5$', fontsize=16)
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(loc='upper right', fontsize=16)
    plt.savefig('convergence.pdf', bbox_inches = 'tight')
    plt.show()

def l2err(u, uex):
    err = 0
    for i in range(u.shape[0]):
        err += (u[i] - uex[i])**2
    err = np.sqrt(1/u.shape[0] * err)

    return err

def characteristics(region):
    ts = np.array([0.5, 1.0, 1.5]); ints = 10

    if region == '1':
        x = np.linspace(0, 1, num=ints, endpoint=True)
        for i in range(x.shape[0]):
            plt.plot([x[i],x[i]],[ts[0], ts[2]], lw=2, color='g')
    elif region == '2':
        x = np.linspace(1, 2, num=ints, endpoint=True)
        for i in range(1, x.shape[0]):
            plt.plot([x[i],x[i] + ts[2]*(x[i]-1)],[ts[0], ts[2]], lw=2, color='r')
    elif region == '3':
        x = np.linspace(2, 3, num=ints, endpoint=True)
        for i in range(1, x.shape[0]):
            plt.plot([x[i],x[i] + ts[2]*(3 - x[i])],[ts[0], ts[2]], lw=2, color='m')
    else:
        x = np.linspace(3, 4, num=ints, endpoint=True)
        for i in range(1, x.shape[0]):
            plt.plot([x[i],x[i]],[ts[0], ts[2]], lw=2, color='b')

def analytical(Nx, t):
    x = np.linspace(0, 4, Nx+1, endpoint=True)
    u = np.zeros(Nx + 1)

    x0 = state_init(x, Nx)
    if t != 1.5:
        for i in range(Nx+1):
            if x[i] >= 0 and x[i] <= 1:
                u[i] = 0
            elif x[i] > 1 and (x[i] + t)/(1 + t) - 1 < (3 - (x[i] - 3*t)/(1 - t)):
                u[i] = (x[i] + t)/(1 + t) - 1
            elif x[i] >= 2 and x[i] < 3:
                if (1 - t) == 0:
                    u[i] == 0
                else:
                    u[i] = 3 - (x[i] - 3*t)/(1 - t)
            elif x[i] >= 3 and x[i] <=4:
                u[i] = 0
    else:
        for i in range(Nx+1):
            if x[i] >= 0 and x[i] <= 1:
                u[i] = 0
            elif x[i] > 1 and x[i] <= 1+np.sqrt(2)*np.sqrt(1+t):
                u[i] = (x[i] + t)/(1 + t) - 1
            else:
                u[i] = 0
    return x, u

def state_init(x, Nx):
    x0 = np.zeros(Nx + 1)
    for i in range(Nx+1):
        if x[i] == 0 or x[i] < 1:
            x0[i] = 0
        elif x[i] == 1 or x[i] < 2:
            x0[i] = x[i] - 1
        elif x[i] == 2 or x[i] < 3:
            x0[i] = 3 - x[i]
        else:
            x0[i] = 0

    return x0

def getIC(Nx):
    x = np.linspace(0, 4, Nx + 1)
    xc = 0.5*(x[0:Nx] + x[1:Nx+1])

    u = np.zeros(Nx+1)
    for i in range(Nx+1):
        if x[i] == 0 or x[i] < 1:
            u[i] = 0
        elif x[i] == 1 or x[i] < 2:
            u[i] = xc[i] - 1
        elif x[i] == 2 or x[i] < 3:
            u[i] = 3 - xc[i]
        else:
            u[i] = 0

    return x, u

if __name__ == "__main__":
    q1a()
    q1b()
    q1c()