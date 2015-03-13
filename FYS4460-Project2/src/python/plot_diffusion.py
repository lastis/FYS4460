from scitools.std import *
import matplotlib.pyplot as plt

filename = '../../res/Measurements/Diffusion.dat'
inFile = open(filename, "r");

lines = inFile.readlines();
n = len(lines);

r2 = zeros(n)
for i in range(n):
    r2[i] = float(lines[i]);
dumpRate = 100;
dt = 0.2;
t = array(range(n))*dt*dumpRate;


plt.figure()
plt.plot(t, r2);
plt.show();

plt.figure()
plt.plot(t[1:-1], r2[1:-1]/(6*t[1:-1]));
plt.show();