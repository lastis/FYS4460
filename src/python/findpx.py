from pylab import*
from scipy import interpolate
import numpy as np

def findNearest(array, value):
	i = (abs(array-value)).argmin()
	return array[i]

files 	= ["25.npz", "50.npz","100.npz","150.npz","200.npz","300.npz",\
		"400.npz","500.npz","800.npz"]
l	= [25.,50.,100.,150.,200.,300.,400.,500.,800.]
px1 	= zeros(len(files))
px2 	= zeros(len(files))
x1 	= 0.8
x2	= 0.2

i = 0
for s in files:
	npzfile = np.load(s)
	p = npzfile['arr_0']
	Pi = npzfile['arr_1']
	
	tck 	= interpolate.splrep(p,Pi,k=3)
	pNew 	= linspace(p[0],p[-1],1001)
	PiNew 	= interpolate.splev(pNew,tck)


	px1[i] = pNew[where(PiNew == findNearest(PiNew,x1))]
	px2[i] = pNew[where(PiNew == findNearest(PiNew,x2))]

	i += 1
	plot(pNew,PiNew)
show()

plot(log(l),log(px1-px2),"-o")
show()


