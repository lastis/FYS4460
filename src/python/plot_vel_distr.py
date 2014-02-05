from scitools.std import *
import matplotlib.pyplot as plt

filename = "/../../res/result/000.xyz"
inFile = open(filename, 'r');
lines = inFile.readlines();
inFile.close();
nlines = len(lines);

#Arrays for velocities u=xvel v=yvel w=zvel
u = zeros(nlines);
v = zeros(nlines);
w = zeros(nlines);

for i in range(2,nlines):
    parts = lines[i].split();
    u[i] = float(parts[4]);
    v[i] = float(parts[5]);
    w[i] = float(parts[6]);
    
# Check mean and standard deviation
print ("x-velocity: mean=%f std=%f" % (u.mean(), u.std()))
print ("y-velocity: mean=%f std=%f" % (v.mean(), v.std()))
print ("z-velocity: mean=%f std=%f" % (w.mean(), w.std()))


def velHistPlot(samples, title):
    plt.figure()
    plt.hist(samples, 40);
    plt.xlabel('Velocity [Aangstrom/s]')
    plt.ylabel('Number of atoms with given velocity')
    plt.title(title);

velHistPlot(u, 'x-velocity distribution')
velHistPlot(v, 'y-velocity distribution')
velHistPlot(w, 'z-velocity distribution')
plt.show()
raw_input()
