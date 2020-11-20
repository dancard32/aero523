import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri, writegri
from plotmesh import plotmesh

def main ():
    mesh = readgri('test1.gri')
    V = mesh['V']; E = mesh['E']; BE = mesh['BE']; IE = mesh['IE']

    print(V, E, BE, IE)
    
    plotmesh(mesh, 'test')

if __name__ == "__main__":
    main()