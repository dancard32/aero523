import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri

#-----------------------------------------------------------
def plotmesh(Mesh, fname):
    V = Mesh['V']; E = Mesh['E']; BE = Mesh['BE']

    f = plt.figure(figsize=(12,12))
    #plt.tripcolor(V[:,0], V[:,1], triangles=E)
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='black')
    dosave = not not fname
    plt.axis('equal')
    plt.axis('off')
    plt.tick_params(axis='both', labelsize=12)
    f.tight_layout(); plt.show(block=(not dosave))
    if (dosave): plt.savefig(fname, bbox_inches='tight')
    plt.close(f)
    
#-----------------------------------------------------------
def main():
    Mesh = readgri('mesh0.gri')
    plotmesh(Mesh, [])

if __name__ == "__main__":
    main()
