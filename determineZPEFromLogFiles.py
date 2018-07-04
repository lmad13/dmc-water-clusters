import numpy as np
import sys
filename=sys.argv[1]
logfile=open(filename,'r')
E=[]

for line in logfile:
    if "average from this simulation" in line:
        E.append(float(line.split()[4]))

E=np.array(E)

print 'E:',E

print ''

print '<E>',np.average(E), 'standard deviation', np.std(E)
