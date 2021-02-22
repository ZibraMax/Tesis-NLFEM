import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath #Para pip

# Adaptado de https://www.fundza.com/vectors/point2line/index.html
def pnt2line(punto, inicio_linea, fin_linea):
	line_vec = fin_linea - inicio_linea
	pnt_vec =  punto - inicio_linea
	line_len = np.linalg.norm(line_vec)
	line_unitvec = line_vec/line_len
	pnt_vec_scaled = pnt_vec * 1.0/line_len
	t = np.dot(line_unitvec, pnt_vec_scaled)  
	if t < 0.0:
		t = 0.0
	elif t > 1.0:
		t = 1.0
	nearest = line_vec * t
	dist = np.linalg.norm(pnt_vec-nearest)
	return dist

def A(g,elemento,xnl,ynl):
	polygon = mpltPath.Path(elemento)
	pip = polygon.contains_points([[xnl,ynl]])
	if not pip[0]:
		d = []
		for i in range(len(coords)-1):
			punto = np.array([xnl,ynl])
			inicio_linea = np.array([elemento[i][0],elemento[i][1]])
			fin_linea = np.array([elemento[i+1][0],elemento[i+1][1]])
			d.append(pnt2line(punto, inicio_linea, fin_linea))
		punto = np.array([xnl,ynl])
		inicio_linea = np.array([elemento[-1][0],elemento[-1][1]])
		fin_linea = np.array([elemento[0][0],elemento[0][1]])
		d.append(pnt2line(punto, inicio_linea, fin_linea))
		return g(np.min(d))
	return g(0)


l = 0.1 #Longitud interna
LR = 6*l
macauley = lambda f, *args: f(*args)/2 + abs(f(*args))/2
f2p = lambda r: (1-(r/LR))
f3p = lambda r: (1-(r**2/LR**2))
g1 = lambda r: np.exp(-r/l) #Función de atenuación biexponencial
g2 = lambda r: macauley(f2p,r)
g3 = lambda r: macauley(f3p,r)
g4 = lambda r: r*g1(r)

G = [g1,g2,g3,g4]
# r_clasico = lambda xl,yl,xnl,ynl: np.sqrt((xnl-xl)**2+(ynl-yl)**2)
r = 1 # Radio del polígono
N = 150 #Numero de puntos para graficar
n = 6 #Numero de lados del poligono
h = 2*(r+LR)/N
th = 2*np.pi/n
coords = np.array([[r*np.cos(th*i),r*np.sin(th*i)] for i in range(n)]) #Coordenadas del polígono
_coords = np.array(coords.tolist()+[coords[0].tolist()]) #Coordenadas del polígono para graficar
_X = [-(r+LR)+h*i for i in range(N+1)] # Muestreo en X
_Y = [-(r+LR)+h*i for i in range(N+1)] # Muestreo en Y
fig = plt.figure()
count = 0
for i,g in enumerate(G):
	count+=1
	X = []
	Y = []
	Z = []
	for x in _X:
		for y in _Y:
			X.append(x) # Coordenada x
			Y.append(y) # Coordenada y
			z = A(g,coords,x,y) # Valor de la función de atenuación
			Z.append(z)
	ax = fig.add_subplot(int(np.sqrt(len(G))),int(np.sqrt(len(G))),count,projection='3d')
	surf = ax.plot_trisurf(X,Y,Z,cmap='magma') # Gráfica 3D
	cbar = fig.colorbar(surf)

plt.show()
