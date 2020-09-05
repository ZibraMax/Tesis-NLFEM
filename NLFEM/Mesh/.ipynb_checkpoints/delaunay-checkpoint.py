import triangle as tr

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss as roots

from mpl_toolkits import mplot3d
from matplotlib import tri
import matplotlib as mpl
from matplotlib.lines import Line2D
from ..Utils import isBetween
import pandas as pd

class Geometria:
    def __init__(this, vertices):
        this.vertices = vertices
        this.cbe = []
        this.cbn = []
class Delaunay1V(Geometria):
    def __init__(this, vertices, params,cbe = [],valorCBE=10, plot=False):
        this.bordesCBE = cbe
        this.valorCBE= valorCBE 
        this.params = params
        super().__init__(vertices)
        this.seg = []
        for i in range(len(this.vertices)-1):
            this.seg.append([i,i+1])
        this.seg.append([i+1,0])
        this.original = dict(vertices=np.array(this.vertices),segments=np.array(this.seg))
        this.triangular = tr.triangulate(this.original,this.params)
        if plot:
            tr.compare(plt, this.original, this.triangular)
        count = 0
        for i in this.triangular['segments']:
            if count > 1:
                if np.sum(np.isin(np.array(this.cbe)[:,0], i[0]))<1:
                    this.cbe.append([i[0],0])
                if np.sum(np.isin(np.array(this.cbe)[:,0], i[1]))<1:
                    this.cbe.append([i[1],0])
            else:
                this.cbe.append([i[0],0])
            if np.sum(np.isin(np.array(this.cbe)[:,0], i[1]))<1:
                    this.cbe.append([i[1],0])
            count+=1
        this.diccionarios = this.triangular['triangles'].tolist()
        this.tipos = np.zeros([len(this.diccionarios)]).astype(str)
        this.tipos[:] = 'T1V'
        this.gdls = this.triangular['vertices'].tolist()
    def darNodosCB(this,segmento):
        a = []
        ps=this.original['vertices'][this.original['segments'][segmento]].tolist()
        for i,p in enumerate(this.triangular['vertices']):
            if isBetween(ps[0], ps[1], p):
                a.append(i)
        return np.array(a)
    def generarCB(this,bordes,valor=0):
        cb = []
        this.bordesCBE = bordes
        this.valorCBE= valor
        for segmento in bordes:
            nodos = this.darNodosCB(segmento)
            cbe = np.zeros([len(nodos),2])
            cbe[:,0] = nodos
            cbe[:,1] = valor
            cb+= cbe.tolist()
        return cb
    def generarDatos(this):
        m = []
        for t in this.triangular['triangles'].tolist():
            coordsModelo = this.original['vertices'].flatten().tolist()
            # Coordendas de cada triangulo
            coords = np.array(this.gdls)[np.ix_(t)]
            cx = [np.average(coords[:,0])]
            cy = [np.average(coords[:,1])]
            coords = coords.flatten().tolist()
            #Saber si es cb
            cb = [np.any(np.isin(t,np.array(this.cbe)[:,0]))*1]
            fila = coordsModelo+coords+cx+cy+cb
            m.append(fila)
        return m
    def areaRefiner(this,model,norm):
        X1 = 'X1'
        X2 = 'Y1'
        X3 = 'X2'
        X4 = 'Y2'
        X5 = 'X2'
        X6 = 'Y3'
        X7 = 'X4'
        X8 = 'Y4'
        X9 = 'X5'
        X10 = 'Y5'
        X11 = 'X6'
        X12 = 'Y6'
        X13 = 'CX'
        X14 = 'CY'
        X15 = 'XE1'
        X16 = 'YE1'
        X17 = 'XE2'
        X18 = 'YE2'
        X19 = 'XE3'
        X20 = 'YE3'
        X21 = 'CB'
        dt = pd.DataFrame.from_records(this.generarDatos())
        dt.columns = [X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X21]
        xs = norm(dt)
        ys = model(xs.values)
        this.triangular['triangle_max_area'] = ys
        tnueva = tr.triangulate(this.triangular,'ra')
        tr.compare(plt,this.triangular,tnueva)
        this.triangular = tnueva
        this.diccionarios = this.triangular['triangles'].tolist()
        this.tipos = np.zeros([len(this.diccionarios)]).astype(str)
        this.tipos[:] = 'T1V'
        this.gdls = this.triangular['vertices'].tolist()
        this.cbe = this.generarCB(this.bordesCBE,this.valorCBE)
        
def Imesh(tf,tw,a,b,params,plot=True):
    corners = [[0,0],[a,0],[a,tf],[a/2+tw/2,tf],[a/2+tw/2,tf+b],[a,tf+b],[a,2*tf+b],[0,2*tf+b],[0,tf+b],[a/2-tw/2,tf+b],[a/2-tw/2,tf],[0,tf]]
    seg = []
    for i in range(len(corners)-1):
        seg.append([i,i+1])
    seg.append([i+1,0])
    A = dict(vertices=np.array(corners),segments=np.array(seg))
    B = tr.triangulate(A,params)
    if plot:
        tr.compare(plt, A, B)
    return B,corners

def _strdelaunay(constrained=True,delaunay=True,a=None,q=None):
    p = ''
    if constrained:
        p = 'p'
    if a == None:
        a = ''
    else:
        a = 'a'+format(a)
    D = ''
    if delaunay:
        D = 'D'
    if q == None:
        q=''
    else:
        if type(q) == int:
            if q > 35:
                raise "No sepuedecrearunatriangulacion conangulos menores a 35 grados"
        q = 'q'+format(q)
    return p+a+D+q+'i'
def generarGeometria(triang):
    cbe = []
    count = 0
    for i in triang['segments']:
        if count > 1:
            if np.sum(np.isin(np.array(cbe)[:,0], i[0]))<1:
                cbe.append([i[0],0])
            if np.sum(np.isin(np.array(cbe)[:,0], i[1]))<1:
                cbe.append([i[1],0])
        else:
            cbe.append([i[0],0])
        if np.sum(np.isin(np.array(cbe)[:,0], i[1]))<1:
                cbe.append([i[1],0])
        count+=1
    g = Geometria([-1])
    
    g.triangular = triang
    g.diccionarios = triang['triangles'].tolist()
    g.gdls = triang['vertices'].tolist()
    g.tipos = np.zeros([len(g.diccionarios)]).astype(str)
    g.tipos[:] = 'T1V'
        
    return g
    