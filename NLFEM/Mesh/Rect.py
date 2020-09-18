import numpy as np
from .Geometria import *

class Rect(Geometria):
	def __init__(this,file,nx=30,ny=30):
		gdls = []
		diccionarios = []
		diccionariosnl = []
		vertices = [-1]
		segmentos = []
		with open(file) as archivo:
			params = archivo.readline().split('	 ')
			ngdls = int(params[0])
			nele = int(params[1])
			tipos = np.zeros([nele]).astype(str)
			tipos[:] = 'C1V'
			for _ in range(ngdls):
				linea = archivo.readline().split('	   ')
				gdls.append([float(linea[0]),float(linea[1])])
			for _ in range(nele):
				linea = list(map(lambda x: int(x)-1,archivo.readline().split('	 ')))
				diccionarios.append(linea)
			for _ in range(nele):
				linea = list(map(lambda x: int(x)-1,archivo.readline().split('	')))
				diccionariosnl.append(linea)
		segmentos.append([0,diccionarios[nx-1][1]])
		segmentos.append([diccionarios[nx-1][1],diccionarios[nx*ny-1][2]])
		segmentos.append([diccionarios[nx*ny-1][2],diccionarios[nx*ny-nx-1][3]])
		segmentos.append([diccionarios[nx*ny-nx-1][3],0])
		super().__init__(vertices, diccionarios, gdls, tipos, segmentos=segmentos)
		this.diccionariosnl = diccionariosnl