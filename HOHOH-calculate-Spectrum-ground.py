#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import os
import glob
import usefulFunctions as use


au2wn=219474.63
nBins=51

if len(sys.argv)<4:
    print 'Usage: ./HOHOH-groundState.py N_size nReps descendantSteps nRepsDW'
    end


starttime=time.time()

stateGround='stateGround'
state='stateAsymSt'
DWstate='DWAsymSt'
molecule='H3O2'
dTau=10

N_size=int(sys.argv[1])
nReps=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])

#figure out which nrep number we're at in the directory of interest
fileParameterName=molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)
path='data'+molecule+'/'+state+'/'+DWstate+'/'
print 'path: ', path
preExistingFiles=glob.glob(path+'*'+fileParameterName+'*')
print 'the files that already exist are:', preExistingFiles

outputFile=open(path+fileParameterName+'-logFile.data','w')
outputFile.write('the files that already exist are: '+str(preExistingFiles)+'\n')

equilibrationSteps=500
propagationSteps=500
averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size)
Wfn.setNodalSurface('OHStretchAnti','Both')

#for each iwfn in nReps, 
for iwfn in range(nReps):
    print '   REPETITION NUMBER: ', iwfn
    groundPath='data'+molecule+'/stateGround/DWGround/'

    #load in the ground state
    groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    GroundCoords,groundDW=Wfn.loadCoords(groundPath+groundStateWfnName)

    eckRotcoords=Wfn.molecule.eckartRotate(GroundCoords)
    nwalkers=GroundCoords.shape[0]

    symCoords,symDW=Wfn.molecule.symmetrizeCoordinates(GroundCoords,groundDW)

    symEckRotCoords=Wfn.molecule.eckartRotate(symCoords)

    GfileName='TheGMatrix-symmetrized'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
    
    #GfileName='TheGMatrix.gmat'
    
    print 'INTERNAL COORDINATES :-O'
    internals=Wfn.molecule.SymInternals(symEckRotCoords)


    print zip(Wfn.molecule.internalName,np.average(internals,weights=symDW,axis=0)*Wfn.molecule.internalConversion)
    
    Wfn.LoadG(groundPath+GfileName,symEckRotCoords,symDW)


    Wfn.calculateSpectrum(symEckRotCoords,symDW)
