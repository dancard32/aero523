import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri

#-----------------------------------------------------------
def plotmesh(Mesh, fname):
    V = Mesh['V']; E = Mesh['E']; BE = Mesh['BE']

    f = plt.figure(figsize=(12,12))
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    for i in range(BE.shape[0]):
        plt.plot(V[BE[i,0:2],0],V[BE[i,0:2],1], '-', linewidth=2, color='black')
    plt.axis('equal'); plt.axis('off')
    f.tight_layout(); 
    plt.savefig(fname, bbox_inches='tight')
    plt.close()
