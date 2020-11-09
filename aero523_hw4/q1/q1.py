import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def main():
    
    plt.figure(figsize=(6,6))
    for t in np.array([0, 0.5, 1.0, 1.5]):
        x0, u0 = analytical(4, t)
        plot_label = r't=' + str(t)
        plt.plot(x0, u0, lw = 2, label=plot_label)
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'Analytical state, $u(x,t)$', fontsize=16)
    plt.legend(loc='upper right', fontsize=16)
    #plt.savefig('analytical.pdf', bbox_inches='tight')
    plt.show()
    

    plt.figure(figsize=(8,6))
    x, slopes = characteristics(4)
    plt.plot(np.NaN, np.NaN, lw=2, color='g')
    plt.plot(np.NaN, np.NaN, lw=2, color='k')
    plt.plot(np.NaN, np.NaN, lw=2, color='m')
    plt.plot(np.NaN, np.NaN, lw=2, color='b')
    plt.plot([0, 0],[-1,1.5], [0.25, 0.25],[-1,1.5], [0.5, 0.5],[-1,1.5], [0.75, 0.75],[-1,1.5], [1,1],[-1,1.5], lw=2, color='g')
    plt.plot([0, 2.5],[-1,1.5], [-2/3, 8/3],[-1,1.5], [-2, 3],[-1,1.5], [-6, 4],[-1,1.5], lw=2, color='k')
    plt.plot([0, 4],[-1, 3], [-2/3,14/3],[-1,3], [-2,6],[-1,3], [-6,4],[-1,1.5], linestyle=':', lw=2, color='m')
    plt.plot([3, 3],[1,3], [3.25, 3.25],[1,3], [3.5, 3.5],[1,3], [3.75, 3.75],[1,3], [4,4],[1,3], lw=2, color='b')
    plt.xlabel(r'x-axis', fontsize=16)
    plt.ylabel(r'Analytical state, $u(x,t)$', fontsize=16)
    plt.xlim((0,4))
    plt.legend([r'Region 1', r'Region 2', r'Region 3', r'Region 4'], loc='lower center', fontsize=16)
    #plt.savefig('analytical.pdf', bbox_inches='tight')
    plt.show()
    
    
def characteristics(Nx):
    x = np.linspace(0, 4, Nx+1, endpoint=True)
    slopes = np.zeros((4, Nx+1))

    slopes[0,:] = 0
    slopes[1,:] = x - 1
    slopes[2,:] = 3 - x
    slopes[3,:] = 0

    return x, slopes

def analytical(Nx, t):
    x = np.linspace(0, 4, Nx+1, endpoint=True)
    u = np.zeros(Nx + 1)

    x0 = state_init(x, Nx)
    for i in range(Nx+1):
        if x[i] == 0 or x[i] < 1:
            u[i] = 0
        elif x[i] == 1 or x[i] < 2:
            u[i] = x0[i] + t*(x0[i] - 1)
        elif x[i] == 2 or x[i] < 3:
            u[i] = x0[i] + t*(3 - x0[i])
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

if __name__ == "__main__":
    main()