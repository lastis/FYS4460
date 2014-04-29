from pylab import*
from scipy.ndimage import measurements
import numpy as np


x 	= 0.8
pc 	= 0.59275
nsample = 5000
p 	= arange(0.45,0.8,0.01)
l      	= array([150.,300.,500.])
ll      = l*l

Pi = zeros((len(l),len(p)))

for j in xrange(0,len(l)):
        sys.stdout.write("\r%s\r" %(" "*30))
        print "Generating for L = " + str(l[j])
        for sample in xrange(0,nsample):
                for i in xrange(0,len(p)):
                        # Print to terminal every 1 percen
                        if(sample%int(nsample/100) == 0):
                            cent = int(sample*100/nsample)
                            cent20 = cent/5
                            sys.stdout.write("\r")
                            sys.stdout.write("[%-20s] %d%%" %("="*cent20,cent))
                            sys.stdout.flush()
                        #Calculate
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
sys.stdout.write("\r%s\r" %(" "*30))
print "Done" 
Pi /= nsample

for i in xrange(len(l)):
        np.savez(str(int(l[i])), p, Pi[i])
	plot(p,Pi[i])
show()
