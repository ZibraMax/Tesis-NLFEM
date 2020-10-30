import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import Rect
import os
PATH = os.getcwd()
nombre = PATH+"\\NLFEM\\Mesh\\input.txt"
GEOMETRIA = Rect.Rect(nombre,10,10)
GEOMETRIA.segmentos = [[0,60],[60,2820],[2820,2760],[2760,0]]
GEOMETRIA.cbe = []
nodos = [0,1,2,3,4,5]
print(np.array(GEOMETRIA.gdls)[np.ix_(nodos)][:,1])

