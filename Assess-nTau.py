#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import os
import glob
import usefulFunctions as use

au2wn=219474.63
nBins=51
if len(sys.argv)<6:
    print 'Usage: ./thisScript.py MOLECULE={H3O2,D3O2} STATE={AsymSt,DeltaR,Ground,ZDisp} N_size nReps descendantSteps nRepsDW '
    end


stateGround='stateGround'

molecule=sys.argv[1]
state='state'+sys.argv[2]
DWstate='DW'+sys.argv[2]
N_size=int(sys.argv[3])
nReps=int(sys.argv[4])
descendantSteps=int(sys.argv[5])
nRepsDW=int(sys.argv[6])

dTau=10
Wfn=dmc.wavefunction('HOHOH', N_size)

excitedPath='data'+molecule+'/'+state+'/'+DWstate+'/'
print 'path: ', excitedPath

groundPath='data'+molecule+'/stateGround/DWGround/'
print 'Ground Path: ', groundPath

groundToExcPath='data'+molecule+'/'+stateGround+'/'+DWstate+'/'
print 'Ground to Excited State path',groundToExcPath

Wfn=dmc.wavefunction('HOHOH', N_size)
if 'D' in molecule:
    Wfn.setIsotope('fullyDeuterated')
    print 'setting deuteration to fully deuterated!'

if sys.argv[2]=='DeltaR':
    surfaceName='SharedProton'
elif sys.argv[2]=='AsymSt':
    signFunc=Wfn.molecule.calcStretchAnti
    surfaceName='OHStretchAnti'
elif sys.argv[2]=='ZDisp':
    surfaceName='Z-displacement'
else:
    print '\nINVALID CHOICE of STATE!\n'
    print 'Choose from: {AsymSt,DeltaR,Ground,ZDisp}'
    die


starttime=time.time()

#figure out which nrep number we're at in the directory of interest
fileParameterName=molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)
path='data'+molecule+'/'+state+'/'+DWstate+'/'
print 'path: ', path
preExistingFiles=glob.glob(path+'*'+fileParameterName+'*')
print 'the files that already exist are:', preExistingFiles

Wfn=dmc.wavefunction('HOHOH', N_size)
Wfn.setNodalSurface(surfaceName,'Both')

#for each iwfn in nReps, 
for iwfn in range(nReps):
    print '   REPETITION NUMBER: ', iwfn
    groundPath='data'+molecule+'/stateGround/DWGround/'

    #load in the ground state
    groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    GroundCoords,groundDW=Wfn.loadCoords(groundPath+groundStateWfnName)


    #propagate for DW  G2E=GroundToExcited
    descendantWeightsG2E=np.zeros((GroundCoords.shape[0]))
    for ides in range(nRepsDW):
        print 'DW Rep Number',ides,
        v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(GroundCoords,descendantSteps,initialPop=N_size,printCensus=False)
        descendantWeightsG2E=descendantWeightsG2E+descendantsTemp
        np.savetxt(path+'vref-pop-'+str(iwfn)+'-'+fileParameterName+'-ides-'+str(ides)+'-array.data',np.array(zip(np.arange(descendantSteps+1),v_ref_DW_list,pop_DW_list)))
        print Wfn.calculateInverseParticipationRatio(descendantsTemp)
    print ''
    descendantWeightsG2E=descendantWeightsG2E/nRepsDW

    print Wfn.calculateInverseParticipationRatio(descendantWeightsG2E)


print 'done! That took',time.time()-starttime,'s'
