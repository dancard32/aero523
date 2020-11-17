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
    u0 = getIC(alpha, 2.2)
    n = np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])

    # Consistency Check
    F, analytical, ignore = RoeFlux(u0, u0, n, True); diff = abs(F - analytical)
    print('Roe Flux Tests:\nConsistency Check\n' + 50*'-' + '\n', F,'\n', analytical)

    f = open('q1/consistency', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$F(u_l,u_l,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(F[0],F[1],F[2],F[3]))
    f.write(r'$\vec{F}(\vec{U}_l)\cdot \vec{n} $ & %.3f & %.3f & %.3f & %.3f \\'%(analytical[0],analytical[1],analytical[2],analytical[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Flipping with Direction
    ul = getIC(alpha, 2.2); ur = getIC(alpha, 2.4); 
    Fl = RoeFlux(ul,ur, n, False); Fr = -1.0*RoeFlux(ur,ul, -n, False); diff = abs(Fl-Fr)
    print('\n\nFlipping Direction\n' + 50*'-' + '\n', Fl,'\n', Fr)

    f = open('q1/flipped', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$F(u_l,u_r,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fl[0],Fl[1],Fl[2],Fl[3]))
    f.write(r'$-F(u_r,u_l,-\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fr[0],Fr[1],Fr[2],Fr[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Supersonic Normal Velocity
    Fl, lanalytic, ranalytic = RoeFlux(ul,ur, n, True)
    Fl1, lanalytic1, ranalytic1 = RoeFlux(ul,ul, n, True)
    Fl2, lanalytic2, ranalytic2 = RoeFlux(ur,ur, n, True)


    
    
def run_fvm():

    solve(0)

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
    test_flux()

    """
    run_fvm()
    mesh_adapt()
    vary_alpha()
    """