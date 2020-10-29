import NLFEM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath

e = NLFEM.CuadrilateroL([[-1.5,0],[1,0],[2,2],[0,1]],[0,1,2,3])
U = np.array([[10],[12],[13],[9]])
x = 0.25
y = 0.25
e.darSolucion(U, graficar=True)
print(e.estaDentro(x,y),e.darSolucionXY( U, x, y, n=100))
