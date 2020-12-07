from .FEM1V import FEM1V
from .. import Mesh
import numpy as np
import matplotlib.pyplot as plt
import triangle as tri
from .polygonal import generatePolygon

def L1error(modeloL,modeloH,plot=False):
    Ug = modeloL.U
    Uf = []
    for gdl in modeloL.geometria.gdls:
        u = None
        for e in modeloH.elementos:
            u = e.solucionInterpoladaLocal(gdl[0],gdl[1])
            if u:
                break
        if not u:
            raise Exception(gdl)
            u = 10
        Uf.append(u)
    Uf = np.array(Uf).reshape([len(Uf),1])
    error = lambda u,u_gt: np.abs((u-u_gt)/u_gt)
    errores = error(Ug,Uf)
    if plot:
        modeloL.graficarSolucionFast(errores,figsize=None,linewidth=0.1,markersize=0.1,name='Distribucion del error')
    return errores,Uf,Ug

def problema(area,points,bcb=[0,1],value=10,plot=False):
    G = 1
    params = Mesh.delaunay._strdelaunay(constrained=True,delaunay=True,a=area)
    vertices = points
    geometria = Mesh.Delaunay1V(vertices, params)
    zanahorias = FEM1V(geometria)
    zanahorias.generarElementos()

    a11 = lambda x,y: 1
    a12 = lambda x,y: 0
    a21 = lambda x,y: 0
    a22 = lambda x,y: 1
    a00 = lambda x,y: 0

    theta = 1
    f = lambda x,y: 2*G*theta
    geometria.cbe = geometria.generarCB(bcb,value)
    zanahorias.definirCondicionesDeBorde(geometria.cbe)
    zanahorias.solucionar(cmap='magma',markersize=1,linewidth=1,a11=a11,a12=a12,a21=a21,a22=a22,a00=a00,f=f,plot=plot)
    return zanahorias

def generarPares(geom=5,radio=5,vert=6,da1=3,da2=100,bcb=[0,1],value=10,plot=False):
    geometrias = []
    for i in range(geom):
        spk = np.random.random()
        irre = np.random.random()
        vertices = generatePolygon( 2*radio, 2*radio, radio, irregularity=spk, spikeyness=irre, numVerts=vert)
        intentos = [radio/da1,radio/da2]
        finitos = []
        for area in intentos:
            problemaOBJ = problema(area,vertices,bcb,value=value,plot=plot)
            if plot:
                problemaOBJ.graficarSolucionFast(problemaOBJ.U,linewidth = 0.1,markersize= 0.1,figsize=None)
                ax = plt.gca()
                ax.set_title('Spikeyness = ' + format(np.round(spk,3)) + ', Irregularity = ' + format(np.round(irre,3)) + ', Area = ' + format(np.round(area,3)))
            finitos.append(problemaOBJ)
        geometrias.append(finitos)
    return geometrias
def areaCorrection(ModeloH,errores,K=0.5, alpha=1,area0=1.67,plot=False):
    
    Areas = K*1/errores**alpha*(K/errores**alpha<=area0)+(K/errores**alpha>area0)*area0
    Areas = np.nan_to_num(Areas, nan=area0)
    if plot:
        ModeloH.graficarSolucionFast(Areas,figsize=None,linewidth=0.1,markersize=0.1,name='Correccion de areas')
    return Areas
def areasRefinado(modeloL,areas):
    a = []
    for e in modeloL.elementos:
        Ue = areas[np.ix_(e.gdl)]
        a.append(np.min(Ue))
        e.areaOptima = a[-1]
    return a