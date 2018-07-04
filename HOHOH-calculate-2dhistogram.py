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
import CalculateSpectrum

au2wn=219474.63
nBins=51

if len(sys.argv)<4:
    print 'Usage: ./HOHOH-groundState.py N_size nReps descendantSteps nRepsDW'
    end


starttime=time.time()

stateGround='stateGround'
state='stateGround'
DWstate='DWGround'
molecule='H3O2'
dTau=10

N_size=int(sys.argv[1])
nReps=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])
nStart=int(sys.argv[5])
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
gatheredSymEckRotCoords=[]
gatheredSymDW=[]
gather2D=np.zeros((51,51))
nWalkersTotal=0
#for each iwfn in nReps, 
for iwfn in range(nStart,nReps):
    print '   REPETITION NUMBER: ', iwfn
    groundPath='data'+molecule+'/stateGround/DWGround/'

    #load in the ground state
    groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    GroundCoords,groundDW=Wfn.loadCoords(groundPath+groundStateWfnName)

    eckRotcoords=Wfn.molecule.eckartRotate(GroundCoords)
    nwalkers=GroundCoords.shape[0]

    symCoords,symDW=Wfn.molecule.symmetrizeCoordinates(GroundCoords,groundDW)

    symEckRotCoords=Wfn.molecule.eckartRotate(symCoords)
    nWalkersTotal+=symCoords.shape[0]

    if iwfn==nStart:
        gatheredSymEckRotCoords=symEckRotCoords
        gatheredSymDW=symDW
    else:
        gatheredSymEckRotCoords=np.append(gatheredSymEckRotCoords,symEckRotCoords,axis=0)
        gatheredSymDW=np.append(gatheredSymDW,symDW,axis=0)

    #print 'two symmeterized walkers!'
    #i=26
    #print symEckRotCoords[i],'\n',symEckRotCoords[i+nwalkers]
    
    print 'INTERNAL COORDINATES :-O'
    internals=Wfn.molecule.SymInternals(symEckRotCoords)
    #print internals[i],'\n',internals[i+nwalkers]
    print np.average(internals,weights=symDW,axis=0)

    outOfPhaseBend=internals[:,4]
    ZDisp=internals[:,8]
    
    histogram2D, xedges,yedges=np.histogram2d(outOfPhaseBend,ZDisp,bins=51, weights=symDW,normed=True)
    
    X_bin_center=(xedges[:-1]+xedges[1:])/2.0
    Y_bin_center=(yedges[:-1]+yedges[1:])/2.0

    gather2D=gather2D+histogram2D


gather2D=gather2D/nReps*1.0

np.savetxt('2DHistrogram.matrix',gather2D)
np.savetxt('2DhistogramAxes.data',zip(X_bin_center,Y_bin_center))

print 'all done!'
