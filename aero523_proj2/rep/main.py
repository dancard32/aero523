import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri
from plotmesh import plotmesh

from flux import RoeFlux
from fvm import solve

def getIC(alpha, mach):
    gam = 1.4
    alpha = np.deg2rad(alpha)
    uinf = np.transpose(np.array([1, mach*np.cos(alpha), mach*np.sin(alpha), 1/(gam*(gam-1)) + mach**2/2]))

    return uinf

def test_flux():
    alpha = 0   
    ul = getIC(alpha, 0.8); ur = getIC(alpha, 0.8)
    n = np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])

    # Consistency Check
    F, analytical, ignore = RoeFlux(ul, ul, n, True); diff = abs(F - analytical)
    print('Roe Flux Tests:\nConsistency Check\n' + 50*'-' + '\n', F,'\n', analytical)

    f = open('q1/consistency', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_l,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(F[0],F[1],F[2],F[3]))
    f.write(r'$\vec{F}(\vec{U}_l)\cdot \vec{n} $ & %.3f & %.3f & %.3f & %.3f \\'%(analytical[0],analytical[1],analytical[2],analytical[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Flipping with Direction
    Fl, fl, fr = RoeFlux(ul,ur, n, True); Fr, fl2, fr2 = RoeFlux(ur,ul, -n, True); Fr *= -1; diff = abs(Fl-Fr)
    print('\n\nFlipping Direction\n' + 50*'-' + '\n', Fl,'\n', Fr)

    f = open('q1/flipped', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_r,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fl[0],Fl[1],Fl[2],Fl[3]))
    f.write(r'$-F(u_r,u_l,-\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fr[0],Fr[1],Fr[2],Fr[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Free-stream 
    Fl, fl, fr = RoeFlux(ul,ul, n, True); Fr, fl2, fr2 = RoeFlux(ur,ur, n, True); diff = abs(Fl-Fr)
    print('\n\nFree Stream Test\n' + 50*'-' + '\n', Fl,'\n', Fr)

    # Supersonic Normal Velocity
    alpha = 0
    ul = getIC(alpha, 2.2); ur = getIC(alpha, 2.5)
    F, FL, FR = RoeFlux(ul,ur, np.array([np.sqrt(2)/2,np.sqrt(2)/2]), True)
    print('\n\nSupersonic Normal Velocity\n'+50*'-'+'\n', F,'\n', FL,'\n', FR)

    f = open('q1/supersonic_normal', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_r,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(F[0],F[1],F[2],F[3]))
    f.write(r'$F_L$ & %.3f & %.3f & %.3f & %.3f \\'%(FL[0],FL[1],FL[2],FL[3]))
    f.write(r'$F_R$ & %.3f & %.3f & %.3f & %.3f'%(FR[0],FR[1],FR[2],FR[3]))
    f.close()
     
def run_fvm():

    u0, u, err, V, E, BE, IE = solve()
    P0 = (1.4 - 1)*(u0[:,3] - 0.5*u0[:,0]*((u0[:,1]/u0[:,0])**2 + (u0[:,2]/u0[:,0])**2))
    P = (1.4 - 1)*(u[:,3] - 0.5*u[:,0]*((u[:,1]/u[:,0])**2 + (u[:,2]/u[:,0])**2))

    """
    f = plt.figure()
    plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=P, cmap='jet', shading='flat')
    plt.axis('equal'); plt.axis('off')
    plt.colorbar()
    plt.show()
    """

    # Part 2
    # Read and process mesh files
    # Iterate until L1 residual norm < 10**-5
    # Calculate ATPR - output to q3 table format
    # Monitor and log residual norm and output convergence - output table format to rep/q2/
    
    # Part 3
    # Also output residual and ATPR to table format to rep/q3/
    # Plot convergence of residaul vs. step iterations
    # Plot ATPR output vs. step iterations
    # Plot two figs. (Mach Number and the Total Pressure)

def mesh_adapt():
    pass

    # Plot sequence of adapted meshes
    # Plot two figs. (Mach Number and the Total Pressure) for the finest mesh
    # Plot ATPR output vs. number of cells in a mesh (last ATPR calculation per iteration)

def vary_alpha():
    pass

    # Vary alpha from 0.5:0.5:3 degrees
    # Run same adaptive iterations for each alpha at least 5
    # Plot ATPR from finest mesh vs. alpha and discuss trend


if __name__ == "__main__":
    #test_flux()
    run_fvm()
    mesh_adapt()
    vary_alpha()


    