import numpy as np
import matplotlib.pyplot as plt

PSIS = [lambda z,n: 1/4*(1-z)*(1-n)*(-1-z-n),lambda z,n: 1/4*(1+z)*(1-n)*(-1+z-n),lambda z,n: 1/4*(1+z)*(1+n)*(-1+z+n),lambda z,n: 1/4*(1-z)*(1+n)*(-1-z+n),lambda z,n: 1/2*(1-z**2)*(1-n),lambda z,n: 1/2*(1+z)*(1-n**2),lambda z,n: 1/2*(1-z**2)*(1+n),lambda z,n: 1/2*(1-z)*(1-n**2)]
DZDPSI = [lambda z,n: -1/4*(n-1)*(2*z+n),lambda z,n: -1/4*(n-1)*(2*z-n),lambda z,n: 1/4*(n+1)*(2*z+n),lambda z,n: 1/4*(n+1)*(2*z-n),lambda z,n: (n-1)*z,lambda z,n: -1/2*(n**2-1),lambda z,n: -(n+1)*z,lambda z,n: 1/2*(n**2-1)]
DNDPSI = [lambda z,n: -1/4*(z-1)*(2*n+z),lambda z,n: 1/4*(z+1)*(2*n-z),lambda z,n: 1/4*(z+1)*(2*n+z),lambda z,n: -1/4*(z-1)*(2*n-z),lambda z,n: 1/2*(z**2-1),lambda z,n: -n*(z+1),lambda z,n: -1/2*(z**2-1),lambda z,n: n*(z-1)]

def enmallado(archivo):
	gdls = []
	elementos = []
	elementosnl = []
	segmentos = []

	with open(file) as archivo:
		params = archivo.readline().split('	 ')
		NGL = int(params[0])
		nele = int(params[1])

		for _ in range(ngdls):
			linea = archivo.readline().split('	   ')
			gdls.append([float(linea[0]),float(linea[1])])

		for _ in range(nele):
			linea = list(map(lambda x: int(x)-1,archivo.readline().split('	 ')))
			elementos.append(linea)

		for _ in range(nele):
			linea = list(map(lambda x: int(x)-1,archivo.readline().split('	')))
			elementosnl.append(linea)

	segmentos.append([0,elementos[nx-1][1]])
	segmentos.append([elementos[nx-1][1],elementos[nx*ny-1][2]])
	segmentos.append([elementos[nx*ny-1][2],elementos[nx*ny-nx][3]])
	segmentos.append([elementos[nx*ny-nx][3],0])

	K = np.zeros([NGL,NGL])
	F = np.zeros([NGL,1])
	Q = np.zeros([NGL,1])

	cbe = []
	cbn = []

def ecuacionesPorElemento():
	pass

def ensamblar():
	pass

def aplicarCondicioneBorde():
	pass

def solucionar():
	pass

def graficarSolucion():
	pass

def integrarElemento2D(p,w):
	pass