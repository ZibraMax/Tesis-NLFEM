import subprocess
E = 2.1*10**6
V = 0.2
t = 0.5

def moduloIntegrador(enmallado,npg,E,V,t,l,RUTA_M,tfa=1):
	subprocess.Popen("./index.exe"+" "+ enmallado +" "+format(int(npg))+" "+format(E)+" "+format(V)+" "+format(t)+" "+format(l)+" "+RUTA_M+" "+format(int(tfa)), shell=False)
ls = [0.1,0.3,0.5,0.8]
ate = [1,2,3,4]

for a in ate:
	for l in ls:
		RUTA_M = 'MATRICES'+'_'+format(a)+'_'+format(l)
		moduloIntegrador("input.txt",3,E,V,t,l,RUTA_M,tfa=a)