#/usr/local/bin/python3
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

data = open(sys.argv[1],"r").readlines()

phi = []


for dat in data:
    if "#" not in dat:
        observables = dat.split()
        phi.append(float(observables[2]))

# the histogram of the data
n, bins, patches = plt.hist(phi, 50, density=1, facecolor='green', alpha=0.75)

plt.xlabel(r'$\langle\psi\rangle$')
plt.ylabel('Probability')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()
