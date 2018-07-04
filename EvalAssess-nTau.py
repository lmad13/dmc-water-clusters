import numpy as np
import glob
import matplotlib.pyplot as plt
import sys
import DMCClusters as dmc


if len(sys.argv)<6:
    print 'Usage: ./thisScript.py MOLECULE={H3O2,D3O2} STATE={AsymSt,DeltaR,Ground,ZDisp} N_size nReps descendantSteps nRepsDW '
    end

molecule=sys.argv[1]
state='state'+sys.argv[2]
DWstate='DW'+sys.argv[2]
N_size=int(sys.argv[3])
nReps=int(sys.argv[4])
descendantSteps=int(sys.argv[5])
nRepsDW=int(sys.argv[6])
dTau=10


#identify files
path='data'+molecule+'/'+state+'/'+DWstate+'/'
fileParameterName=molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)

print 'looking for', path+'vref-pop-*-'+fileParameterName+'-ides-*-array.data'

Wfn=dmc.wavefunction('HOHOH', 20000)

#loadtxt
for iwfn in range(nReps):
    print 'files for ', iwfn
    files=glob.glob(path+'vref-pop-*'+str(iwfn)+'-'+fileParameterName+'-ides-*-array.data')    
    print len(files)
    gatherFirstDerivatives=[]
    for f in files:
        data=np.loadtxt(f)
        time=data[:,0]
        v_ref=data[:,1]
        pop=data[:,2]
        y=pop
        smtht,smthy=dmc.use.takeDerivative(time,y,0,noisy=True,smoothOver=5)
        t1,dydx=dmc.use.takeDerivative(time,y,1,noisy=True)
        gatherFirstDerivatives.append(dydx)
        t2,d2ydx2=dmc.use.takeDerivative(time,y,2,noisy=True)
    
        #plt.plot(time,y)
        avgy=np.average(smthy)
        #plt.plot(smtht,smthy-avgy)
        
        #plt.plot(t1,dydx)
        
        #plt.plot(t2,d2ydx2)

    gatherFirstDerivatives=np.array(gatherFirstDerivatives)

    plt.plot(t1[10:],np.average(gatherFirstDerivatives,axis=0)[10:])
plt.show()

#calc 1st derivative

#calc 2nd derivative
