from NLFEM.Mesh import Rect
from NLFEM.Mesh import delaunay
from NLFEM.Mesh import generador
import matplotlib.pyplot as plt

a = 5
vertices = [[0,0],[a,0],[a,a],[0,a]]
GEOMETRIA = delaunay.Delaunay1V(vertices, delaunay._strdelaunay(constrained=True,delaunay=True,a=0.05,q=30,o=2))
GEOMETRIA.detectarNoLocales(6*0.1)
GEOMETRIA.guardarArchivoEnmallado('NLFEM_C++/input.txt')
# GEOMETRIA.animacionNoLocal()
