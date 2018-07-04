import numpy as np
import glob
import sys

if len(sys.argv)<4:
    print 'args: molecule{H3O2 or D3O2} nWalkers nT nDW'
    end

molecule=sys.argv[1]
nWalkers=sys.argv[2]
nt=sys.argv[3]
nDW=sys.argv[4]
regexp='vref-pop-*-'+molecule+'-state*-DW*-dt10-nWalk'+nWalkers+'-nT'+nt+'-nDW'+nDW+'-array.data'
print regexp

fileNames=glob.glob(regexp)
print 'files:\n',fileNames
E=[]
for fn in fileNames:
    data=np.loadtxt(fn)
    E.append(np.average(data[250:,1]))



conv=219474.63
E=np.array(E)
E=E*conv
print 'E:\n',E

print np.average(E),'+/-',np.std(E),'cm-1'

