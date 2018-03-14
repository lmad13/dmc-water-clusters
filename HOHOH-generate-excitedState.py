#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys


if len(sys.argv)<4:
    print 'Usage: ./HOHOH-groundState.py N_size nReps descendantSteps nRepsDW'
    end

N_size=int(sys.argv[1])
nReps=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])

equilibrationSteps=0
propagationSteps=500

starttime=time.time()
au2wn=219474.63
nBins=51
AvePsi2Hist=np.zeros((nBins))
AvePsi2DipHist=np.zeros((nBins))
AvePsi2R2Hist=np.zeros((nBins))
AvePsi2Dip2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size)
Wfn.setNodalSurface('SharedProton',side='Both')
Destination='ResultsH3O2-Tau/'
GatherExpectationRn=[]
GatherExpectationRn2=[]
GatherExpectationMagMu=[]
GatherExpectationMagMu2=[]


initialx=Wfn.x*1.1

#print 'initial V', Wfn.molecule.V([initialx[0]])*au2wn
#print 'equilibrating for ', equilibrationSteps, 'steps (',equilibrationSteps*Wfn.dtau,' au)'
#v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
inputx=initialx

parameterString=str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)
plotFileName=Destination+'Vref-Pop-histogram-Excited'+parameterString
plt.figure(1)

#Sampling of Psi

for iwfn in range(nReps):
    print '\n   REPETITION NUMBER: ', iwfn
    v_ref_list,pop_list,finalCoords,d=Wfn.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)

    averaged_vref.append(np.average(np.array(v_ref_list)*au2wn))
    list_of_pop_list.append(pop_list)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps),np.array(v_ref_list)*au2wn)

    plt.figure(1)
    plt.subplot(312)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps),np.array(pop_list))

    #Rn=Wfn.molecule.calcSharedProtonDisplacement(finalCoords)
    #Dipole=Wfn.molecule.calcDipole(finalCoords)

    for tauDW in [descendantSteps,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
        descendantWeights=np.zeros((finalCoords.shape[0]))
        for ides in range(nRepsDW):
            print 'DW Rep Number',ides,
            v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(finalCoords,descendantSteps,initialPop=N_size,printCensus=False)
            descendantWeights=descendantWeights+descendantsTemp
        descendantWeights=descendantWeights/nRepsDW
        ##inputx=finalCoords
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        Wfn.exportCoords(finalCoords,'Wfn-HOHOH-Tau/HOHOH-ExcitedRn-'+parameterString+'Eq-'+str(iwfn)+'.xyz',descendantWeights)



endtime=time.time()

plt.savefig(plotFileName+'.png')
plt.clf()


print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'

print 'done!'
