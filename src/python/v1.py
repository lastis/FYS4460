from pylab import*
from scipy.ndimage import measurements

nsample = 5
pc = 0.59275
p = arange(pc-0.1,1.0,0.01)

lStart = 9
lEnd = 9
P = zeros((len(p),lEnd-lStart+1))
Pi = zeros((len(p),lEnd-lStart+1))

for lCnt in xrange(lEnd - lStart + 1): 
	lx = pow(2,lStart + lCnt)
	ly = lx
	ll = lx*ly
	for sample in xrange(0,nsample):
		r = rand(lx,ly)
		for i in xrange(0,len(p)):
			z = r < p[i]
			lw, num = measurements.label(z)
			perc_x = intersect1d(lw[0,:],lw[-1, :])
			perc_y = intersect1d(lw[:,0],lw[ :,-1])
			perc_u = union1d(perc_x,perc_y)
		    	perc = perc_u[where(perc_u > 0)]
			if (len(perc) > 0):
				area = measurements.sum(z,lw,index=perc)
				totSum = sum(area)
				P[i,lCnt] += totSum/ll
				Pi[i,lCnt] += 1

Pi /= nsample
P /=  nsample
plot(p,Pi)
show()

plot(p,P)
show()

P = log(P)
p = log(p-pc)

plot(p,P)
#show()
				
