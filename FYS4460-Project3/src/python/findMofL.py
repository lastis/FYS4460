from pylab import*
from scipy.ndimage import measurements

nsample = 50
pc 	= 0.59275
p 	= pc
k       = arange(4,11)
lx      = pow(2,k)
ly      = lx
ll      = lx*ly
massPerc = zeros(len(k))

logbinsize = 1.15
logbinmax = ll


for i in xrange(0,len(k)):
	for sample in xrange(0,nsample):
		r = rand(lx[i],ly[i])
		z = r < p
		lw, num = measurements.label(z)
		perc_x = intersect1d(lw[0,:],lw[-1, :])
		#perc_y = intersect1d(lw[:,0],lw[ :,-1])
		#perc_u = union1d(perc_x,perc_y)

                # Remove the background if percolating
		perc = perc_x[where(perc_x > 0)]
		if (len(perc) > 0):
			area = measurements.sum(z,lw,index=perc)
			totSum = sum(area)
                        massPerc[i] += totSum

massPerc /= nsample
plot(log(lx),log(massPerc))
show()
