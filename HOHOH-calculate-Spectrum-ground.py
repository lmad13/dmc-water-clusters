#!/usr/bin/python
import numpy as np
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

    fileOut=open('symEckRotCoords.xyz', 'w')
    Wfn.molecule.printCoordsToFile(symEckRotCoords,fileOut)
    fileOut.close()
    #print 'two symmeterized walkers!'
    #i=26
    #print symEckRotCoords[i],'\n',symEckRotCoords[i+nwalkers]
    
    print 'INTERNAL COORDINATES :-O'
    internals=Wfn.molecule.SymInternals(symEckRotCoords)
    #print internals[i],'\n',internals[i+nwalkers]
    print np.average(internals,weights=symDW,axis=0)




    GfileName='TheGMatrix-symmetrized'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
    #GfileName='TheGMatrix-symmetrized-all-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
    
    #GfileName='TheGMatrix.gmat'
    
    HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckRotCoords,symDW,path=path)
    #HOASpectrum.calculateG(symEckRotCoords,symDW)
    fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckRotCoords,symDW,path+GfileName)
    fundamentalFile=open(path+'Fundamentals-'+str(iwfn)+'-from-'+fileParameterName+'.data','w')
    for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
        fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
    fundamentalFile.close()
    combinationFile=open(path+'Combinations-'+str(iwfn)+'-from-'+fileParameterName+'.data','w')
    for i in range(combinationBandEnergies.shape[0]):
        for j in range(i):
            combinationFile.write(str(i)+"  "+str(j)+"       "+str(combinationBandEnergies[i,j])+"    "+str(combinationBandIntensities[i,j])+"\n")
    combinationFile.close()

#gatheredSymEckRotCoords=np.array(gatheredSymEckRotCoords).reshape(nWalkersTotal,Wfn.nAtoms,3)
#gatheredSymDW=np.array(gatheredSymDW).reshape(nWalkersTotal)
HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,gatheredSymEckRotCoords,gatheredSymDW,path=path)
#HOASpectrum.calculateG(gatheredSymEckRotCoords,gatheredSymDW)
GfileName='TheGMatrix-symmetrized-all-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'
fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(gatheredSymEckRotCoords,gatheredSymDW,path+GfileName)
fundamentalFile=open(path+'Fundamentals-from-Wfns-'+fileParameterName+str(nStart)+'-to-'+str(nReps)+'.data','w')
for i,(v,intensity) in enumerate(zip(fundamentalEnergies,fundamentalIntensities)):
    fundamentalFile.write(str(i)+"       "+str(v)+"   "+str(intensity)+"\n")
fundamentalFile.close()
combinationFile=open(path+'Combinations-from-Wfns-'+fileParameterName+str(nStart)+'-to-'+str(nReps)+'.data','w')
for i in range(combinationBandEnergies.shape[0]):
    for j in range(i):
        combinationFile.write(str(i)+"  "+str(j)+"       "+str(combinationBandEnergies[i,j])+"    "+str(combinationBandIntensities[i,j])+"\n")
combinationFile.close()

    #print zip(Wfn.molecule.internalName,np.average(internals,weights=symDW,axis=0)*Wfn.molecule.internalConversion)
    
#    Wfn.LoadG(groundPath+GfileName,symEckRotCoords,symDW)



