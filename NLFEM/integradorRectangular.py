from .Elemento import *
from numpy.polynomial.legendre import leggauss as roots

class integradorRectangular(Elemento):
    def __init__(this,coords,gdl=None,gauss=4):
        super().__init__(coords=coords,gdl=gdl,gauss=gauss)
        DENSIDAD = 10
        XX,YY = np.meshgrid(np.linspace(-1,1,10), np.linspace(-1,1,DENSIDAD))
        this._dominioNaturalZ = XX.reshape([DENSIDAD**2,1])[:,0]
        this._dominioNaturalN = YY.reshape([DENSIDAD**2,1])[:,0]
    def intGauss2D(this,n,f):
        pyp = roots(n)
        X = pyp[0]
        Y = pyp[0]
        Wx = pyp[1]
        Wy = pyp[1]
        INT = 0
        for i,z in enumerate(X):
            for j,n in enumerate(Y):
                INT += f(z,n)*Wx[i]*Wy[j]
        return INT
    def intGauss4DNL(this,n,f,shape=[1]):
        pyp = roots(n)
        X = pyp[0]
        Y = pyp[0]
        Wx = pyp[1]
        Wy = pyp[1]
        INT = np.zeros(shape)
        for i,z in enumerate(X):
            for j,n in enumerate(Y):
                for inl,znl in enumerate(X):
                    for jnl,nnl in enumerate(Y):
                        INT += f(z,n,znl,nnl)*Wx[i]*Wy[j]*Wx[inl]*Wy[jnl]
        return INT
