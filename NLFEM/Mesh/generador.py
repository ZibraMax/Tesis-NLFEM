import numpy as np
import matplotlib.pyplot as plt
# l,lx,ly,nex,ney,filename,animacion = sys.argv[1:]
def meshF(l,lx,ly,nex,ney,filename,animacion):
	l=float(l)
	lx=float(lx)
	ly=float(ly)
	nex=int(nex)
	ney=int(ney)
	lr=6*l
	nne=8
	hx=lx/nex
	hy=ly/ney
	nnd=(ney+1)*(2*nex+1)+(ney)*(nex+1)
	x=np.zeros([nnd])
	y=np.zeros([nnd])
	nel=nex*ney
	nle=[]
	elm=np.zeros([nel,9])
	ntel=1
	#Coordinate Generation
	print('Generando Coordenadas')
	nd=-1
	for i in range(1,ney+1):
		cy=(i-1)*hy
		for j in range(1,2*nex+2):
			nd=nd+1
			y[nd]=cy
			x[nd]=(j-1)*hx/2
		cy=(i-1)*hy+hy/2
		for j in range(1,nex+2):
			nd=nd+1
			y[nd]=cy
			x[nd]=(j-1)*hx
	cy=ly
	for j in range(1,2*nex+2):
		nd=nd+1
		y[nd]=cy
		x[nd]=(j-1)*hx/2
	## Element Node Connectivty    
	print('Generando Elementos')

	ne=-1
	for i in range(0,ney):
		ne=ne+1
		elm[ne,0]=(i)*(3*nex+2)+1
		elm[ne,1]=elm[ne,0]+2
		elm[ne,3]=elm[ne,0]+3*nex+2
		elm[ne,2]=elm[ne,3]+2
		elm[ne,4]=elm[ne,0]+1    
		elm[ne,7]=elm[ne,0]+2*nex+1
		elm[ne,6]=elm[ne,3]+1
		elm[ne,5]=elm[ne,7]+1
		elm[ne,8]=1
		for j in range(1,nex):
			ne=ne+1
			elm[ne,0]=elm[ne-1,1]
			elm[ne,1]=elm[ne,0]+2
			elm[ne,3]=elm[ne-1,2]
			elm[ne,2]=elm[ne,3]+2
			elm[ne,4]=elm[ne,0]+1
			elm[ne,7]=elm[ne-1,5]
			elm[ne,6]=elm[ne,3]+1
			elm[ne,5]=elm[ne,7]+1
			elm[ne,8]=1
	## Identify neighbors within influence distance
	print('Detectando Elementos No Locales')
	nem=0
	for i in range(0,nel):
		ne=-1
		nii=elm[i,0]
		nfi=elm[i,2]
		xci=(x[int(nii-1)]+x[int(nfi-1)])/2
		yci=(y[int(nii-1)]+y[int(nfi-1)])/2  
		ne=ne+1
		g = ne
		nnn = []
		for j in range(0,nel):
			if not j==i:
				nij=elm[j,0]
				nfj=elm[j,2]
				xcj=(x[int(nij-1)]+x[int(nfj-1)])/2
				ycj=(y[int(nij-1)]+y[int(nfj-1)])/2
				dist=np.sqrt((xci-xcj)**2+(yci-ycj)**2)
				if dist<(lr+0.5*np.sqrt(hx**2+hy**2)):
					ne=ne+1
					nnn.append(j+1)
		nle.append(np.array([ne]+[i+1]+nnn))
		if ne>nem:
			nem=ne
	nle = np.array(nle)    
	print('Guardando Archivo')
	f=open(filename,'w')
	f.write(format(nnd)+'\t'+format(nel)+'\n')
	for i in range(nnd):
		f.write(format(x[i])+'\t'+format(y[i])+'\n')
	for i in range(nel):
		fun = lambda x:str(int(x))
		f.write('\t'.join(map(fun,[elm[i,0],elm[i,1],elm[i,2],elm[i,3],elm[i,4],elm[i,5],elm[i,6],elm[i,7]]))+'\n')
	for i in range(nel):
		sring = ''+format(int(nle[i][0]+1))+'\t'
		for j in range(1,len(nle[i])):
			if j==(len(nle[i])-1):
				sep = ''
			else:
				sep = '\t'
			sring+=format(int(nle[i][j]))+sep
		f.write(sring+'\n')
	f.close()
	print('Archivo '+ filename + ' Guardado')
	if not animacion=='0':
		plt.ion()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(x,y,'.',color='black')
		ax.set_aspect('equal')
		for elemento in range(nel):
			for j in range(1,int(nle[elemento][0])+2):
				i = nle[elemento][j]-1
				_x = x[list(map(lambda x: int(x)-1,[elm[i,0],elm[i,1],elm[i,2],elm[i,3]]))]
				_y = y[list(map(lambda x: int(x)-1,[elm[i,0],elm[i,1],elm[i,2],elm[i,3]]))]
				if j==1:
					ax.fill(_x,_y,color='red')
				else:
					ax.fill(_x,_y,color='blue')
			fig.canvas.draw()
			fig.canvas.flush_events()
			ax.patches=[]
