import matplotlib.pyplot as plt
import numpy as np
from scitools.std import *

data = loadtxt('../cpp/fcc.xyz',skiprows=2,usecols=(1,2,3,4,5,6))

r = data[:,[0,1,2]]
v = data[:,[3,4,5]]
vAbs = multiply(v,v)
# Sum the rows
vAbs = vAbs.sum(axis=1)
vAbs = sqrt(vAbs)

hist, bins = np.histogram(vAbs,bins=50)
plt.bar(hist)
plt.show()
