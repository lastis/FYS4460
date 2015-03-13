from pylab import*
from scipy.ndimage import measurements
from numpy import linalg

def logbin(y,a,binmax):
	ymax = binmax
	yedge = 1.0
	istep = 1
	yedgelast = 0
	edge = []
	while (yedgelast <= ymax):
		edge.append(yedge)
		yyedge = floor(yedge*a)
		dy = yyedge - yedge
		if(dy <= 1.0): # We don't want boxes smaller than one
			yyedge += 1
		yedgelast = yedge
		yedge = yyedge
		istep += 1
	n, dummy = histogram(y,edge)
	dx = diff(edge)
	x = 0.5*dx + edge[0:-1]
	return (x,dx,n)


# Variables
nsample = 2
pc 	= 0.59275
p 	= arange(pc-0.2,pc,0.05)


lx = 1500.0 #pow(2,10)
ly = lx
ll = lx*ly

logbinsize = 1.15
logbinmax = 1500

Pi = zeros(len(p))
P = zeros(len(p))
nsp = []

for i in xrange(0,len(p)):
	for sample in xrange(0,nsample):
		r = rand(lx,ly)
		z = r < p[i]
		lw, num = measurements.label(z)
		perc_x = intersect1d(lw[0,:],lw[-1, :])
		perc_y = intersect1d(lw[:,0],lw[ :,-1])
		perc_u = union1d(perc_x,perc_y)
		# This line might be uncesseray
		perc = perc_u[where(perc_u > 0)]
		if (len(perc) > 0):
			Pi[i] += 1
			area = measurements.sum(z,lw,index=perc)
			totSum = sum(area)
			P[i] += totSum/ll
		# Remove percolation cluster
		labels = arange(1,num+1) # might be 0
		delete(labels,perc)
		area = measurements.sum(z, lw, index=labels)
		x, dx, n = logbin(area,logbinsize,logbinmax)

		nnsp = n/ll
		nnsp = nnsp/dx
		if(sample == 0):
			nsp.append(nnsp)
		else :
			nsp[i] += nnsp

Pi /=  nsample
P  /=  nsample

for i in xrange(len(nsp)):
	plot(log(x),log(nsp[i]))
show()

"""
plot(p,P)
show()
"""



	






