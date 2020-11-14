import numpy as np
import matplotlib.pyplot as plt
from readgri import readgri
from plotmesh import plotmesh

def main():
    mesh = readgri('mesh0.gri')
    plotmesh(mesh, 'test.pdf')

if __name__ == "__main__":
    main()