#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
au2wn=219474.63
au2ang=0.529177249
print 'you are going to calculate v for this file', sys.argv[1]

fileName=sys.argv[1]

Wfn=dmc.wavefunction('HOHOH', 100)

finalCoords,descendantWeights=Wfn.loadCoords(fileName)

v=Wfn.molecule.V(finalCoords)
print v
print v*au2wn
print 'if the ZPE of the gorund state is 6604 1/cm',
print '(',6604.0/au2wn,' a.u.)',
print 'and the first excited state is at least 3610 1/cm',
print '(',3610/au2wn,' a.u.)', 
print ', vref would be something like 10,214 1/cm (',(10214.0/au2wn),' a.u.)'
vref=10214/au2wn
for dtau in [1,5,10]:
    print 'rough probabilities for death with dtau=',dtau,' is ',100*(1-np.exp(-(v-vref)*dtau)),'%'


hist,bin_edges=np.histogram(100*(1-np.exp(-(v-vref)*dtau)))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0

plt.plot(bin_center,hist)
plt.show()

