from .integradorTriangular import *
class TriangularL(integradorTriangular):
    def __init__(this,coords,gdl=None,gauss=4):
        if len(coords) == 3:
            this.psis = lambda z,n: np.array([[1-z-n],[z],[n]])
            this.dzpsis = lambda z,n: np.array([[-1],[1],[0]])
            this.dnpsis = lambda z,n: np.array([[-1],[0],[1]])
        else:
            raise Exception('Error: Se esta intentando crear un elemento triangualar que no tiene 3 coordenadas. '+'Recuerde que se necesitan de 3 corrdenadas en orden contrario de las manecillas del reloj '+'sin repetir la primera coordenada')
        super().__init__(coords=coords,gdl=gdl,gauss=gauss)
        corners = np.array(this.coords)
        this.alpha = [corners[1][0]*corners[2][1]-corners[2][0]*corners[1][1],corners[2][0]*corners[0][1]-corners[0][0]*corners[2][1],corners[0][0]*corners[1][1]-corners[1][0]*corners[0][1]]
        this.beta = [corners[1][1]-corners[2][1],corners[2][1]-corners[0][1],corners[0][1]-corners[1][1]]
        this.gamma =[-(corners[1][0]-corners[2][0]),-(corners[2][0]-corners[0][0]),-(corners[0][0]-corners[1][0])]
        this.area2 = np.sum(this.alpha)
        this.psisLocales = lambda x,y: np.array([[1/this.area2*(this.alpha[0]+this.beta[0]*x+this.gamma[0]*y)],[1/this.area2*(this.alpha[1]+this.beta[1]*x+this.gamma[1]*y)],[1/this.area2*(this.alpha[2]+this.beta[2]*x+this.gamma[2]*y)]])
    def solucionInterpoladaLocal(this,x,y):
        t = this.estaDentro(x,y)
        if t:
            return (this.psisLocales(x,y).T@this.Ue)[0,0]
        else:
            return False