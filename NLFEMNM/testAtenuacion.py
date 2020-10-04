from numpy.polynomial.legendre import leggauss as roots

def atenuacion(x,xp,l):
    return 1/(2*l)*np.exp(-np.abs(x-xp)/l)

def intGauss2D(n, f):
    pyp = roots(n)
    X = pyp[0]
    W = pyp[1]
    return sum(f(X, X) @ W * W)
n = 100
l=0.1
xl=-3
xu=3

f = lambda t: intGauss2D(n,lambda x,y : atenuacion(x,y,l,t))-1
for i in range(50):
    xr = (xl + xu)/2
    if f(xr)*f(xl)<0:
        xu = xr
    else:
        xl = xr
    print(i,xr,f(xr))
