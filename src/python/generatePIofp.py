from pylab import*
from scipy.ndimage import measurements

def findNearest(array, value):
	i = (abs(array-value)).argmin()
	return array[i]

x 	= 0.8
pc 	= 0.59275
nsample = 1000
p 	= arange(0.45,0.8,0.01)
l      	= array([25.,50.,100.,200.])
ll      = l*l

Pi = zeros((len(l),len(p)))

for j in xrange(0,len(l)):
	for i in xrange(0,len(p)):
		for sample in xrange(0,nsample):
			r = rand(l[j],l[j])
			z = r < p[i]
			lw, num = measurements.label(z)
			perc_x = intersect1d(lw[0,:],lw[-1, :])
			#perc_y = intersect1d(lw[:,0],lw[ :,-1])
			#perc_u = union1d(perc_x,perc_y)

			# Remove the background if it's percolating
			perc = perc_x[where(perc_x > 0)]
			if (len(perc) > 0):
				Pi[j,i] += 1
Pi /= nsample

for i in xrange(len(l)):
	plot(p,Pi[i])
show()
