import numpy as np
import matplotlib.pyplot as plt
import time

# Project specific functions
from readgri import readgri
from plotmesh import plotmesh
from flux import RoeFlux
from fvm import solve

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
    F, analytical, FR, ls = RoeFlux(ul, ul, n); diff = abs(F - analytical)
    print('Roe Flux Tests:\nConsistency Check\n' + 50*'-' + '\n', F,'\n', analytical)

    f = open('q1/consistency', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_l,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(F[0],F[1],F[2],F[3]))
    f.write(r'$\vec{F}(\vec{U}_l)\cdot \vec{n} $ & %.3f & %.3f & %.3f & %.3f \\'%(analytical[0],analytical[1],analytical[2],analytical[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Flipping with Direction
    Fl, FL, FR, ls = RoeFlux(ul,ur, n); Fr, FL, FR, ls = RoeFlux(ur,ul, -n); Fr *= -1; diff = abs(Fl-Fr)
    print('\n\nFlipping Direction\n' + 50*'-' + '\n', Fl,'\n', Fr)

    f = open('q1/flipped', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_r,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fl[0],Fl[1],Fl[2],Fl[3]))
    f.write(r'$-F(u_r,u_l,-\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(Fr[0],Fr[1],Fr[2],Fr[3]))
    f.write(r'$\Delta F$ & %.2e & %.2e & %.2e & %.2e'%(diff[0],diff[1],diff[2],diff[3]))
    f.close()

    # Free-stream 
    Fl, FL, FR, ls = RoeFlux(ul,ul, n); Fr, FL, FR, ls = RoeFlux(ur,ur, n); diff = abs(Fl-Fr)
    print('\n\nFree Stream Test\n' + 50*'-' + '\n', Fl,'\n', Fr)

    # Supersonic Normal Velocity
    alpha = 0
    ul = getIC(alpha, 2.2); ur = getIC(alpha, 2.5)
    F, FL, FR, ls = RoeFlux(ul,ur, np.array([np.sqrt(2)/2,np.sqrt(2)/2]))
    print('\n\nSupersonic Normal Velocity\n'+50*'-'+'\n', F,'\n', FL,'\n', FR)

    f = open('q1/supersonic_normal', 'w')
    f.write(r'Flux & $\rho$ & $\rho u$ & $\rho v$ & $\rho E$ \\ \hline\hline')
    f.write(r'$\hat{F}(u_l,u_r,\vec{n})$ & %.3f & %.3f & %.3f & %.3f \\'%(F[0],F[1],F[2],F[3]))
    f.write(r'$F_L$ & %.3f & %.3f & %.3f & %.3f \\'%(FL[0],FL[1],FL[2],FL[3]))
    f.write(r'$F_R$ & %.3f & %.3f & %.3f & %.3f'%(FR[0],FR[1],FR[2],FR[3]))
    f.close()
     
def post_process(u):
    gam = 1.4
    uvel = u[:,1]/u[:,0]; vvel = u[:,2]/u[:,0]

    q = np.sqrt(uvel**2 + vvel**2)
    p = (gam-1)*(u[:,3]-0.5*u[:,0]*q**2)
    H = (u[:,3] + p)/u[:,0]    

    c = np.sqrt((gam-1.0)*(H - 0.5*q**2))
    mach = q/c
    pt = p*(1 + 0.5*0.4*mach**2)**(gam/(gam-1))

    return mach, pt

def run_fvm():

    start = time.time()
    u, err, ATPR, V, E, BE, IE = solve(1); end = time.time(); print('Elapsed Time %.2f'%(end - start))
    mach, pt = post_process(u)


    plt.figure(figsize=(8,5))
    plt.plot(np.arange(err.shape[0]), err, lw=2, color='k')
    plt.xlabel(r'Iterations', fontsize=16)
    plt.ylabel(r'$L_1$ Norm', fontsize=16)
    plt.xscale('log'); plt.yscale('log')
    plt.savefig('q3/l1_err.pdf', bbox_inches='tight')
    plt.show()

    plt.figure(figsize=(8,5))
    plt.plot(np.arange(ATPR.shape[0]), ATPR, lw=2, color='k')
    plt.xlabel(r'Iterations', fontsize=16)
    plt.ylabel(r'ATPR Output', fontsize=16)
    plt.savefig('q3/ATPR.pdf', bbox_inches='tight')
    plt.show()

    plt.figure(figsize=(8,4.5))
    plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=mach, vmin=0.9, vmax=2.5, cmap='jet', shading='flat')
    plt.axis('off')
    plt.colorbar(label='Mach Number')
    plt.savefig('q3/Machfield.pdf', bbox_inches='tight')
    plt.show()

    plt.figure(figsize=(8,4.5))
    plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=pt, vmin=6.5, vmax=7.6, cmap='jet', shading='flat')
    plt.axis('off')
    plt.colorbar(label='Total Pressure')
    plt.savefig('q3/Pfield.pdf', bbox_inches='tight')
    plt.show()


def mesh_adapt():

    # Plot sequence of adapted meshes
    # Plot two figs. (Mach Number and the Total Pressure) for the finest mesh
    # Plot ATPR output vs. number of cells in a mesh (last ATPR calculation per iteration)

    mesh = readgri('mesh0.gri')
    for i in range(6):
        plotmesh(mesh, 'q4/mesh' + str(i) + '.pdf')

def vary_alpha():

    # Vary alpha from 0.5:0.5:3 degrees
    # Run same adaptive iterations for each alpha at least 5
    # Plot ATPR from finest mesh vs. alpha and discuss trend
    alphas = np.arange(0.5,3.5, step=0.5)
    atpr_out = np.zeros(6); k = 0
    for i in alphas:

        start = time.time()
        u, err, ATPR, V, E, BE, IE = solve(i); end = time.time(); print('Elapsed Time %.2f'%(end - start))
        mach, pt = post_process(u)

        plt.figure(figsize=(8,4.5))
        plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=mach, vmin=0.9, vmax=2.5, cmap='jet', shading='flat')
        plt.axis('off')
        plt.savefig('q5/mach_a' + str(int(i*10)) + '.pdf', bbox_inches='tight')
        plt.pause(0.2)
        plt.close()

        plt.figure(figsize=(8,4.5))
        plt.tripcolor(V[:,0], V[:,1], triangles=E, facecolors=pt, vmin=6.5, vmax=7.6, cmap='jet', shading='flat')
        plt.axis('off')
        plt.savefig('q5/pt_a' + str(int(i*10)) + '.pdf', bbox_inches='tight')
        plt.pause(0.2)
        plt.close()

        atpr_out[k] = ATPR[len(ATPR)-1]; k += 1


    f = open('q5/atpr_out', 'w'); output = ''
    for i in range(6):
        output += r'%.1f\degree & %.3f \\'%(alphas[i], atpr_out[i])
    f.write(output)
    f.close()

    plt.figure(figsize=(9,5))
    plt.plot(alphas, atpr_out, lw=2, color='k')
    plt.xlabel(r'Angle of attack, $\alpha$', fontsize=16)
    plt.ylabel(r'ATPR Output', fontsize=16)
    plt.savefig('q5/ATPR.pdf', bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    test_flux()
    #run_fvm()
    #mesh_adapt()
    #vary_alpha()
    