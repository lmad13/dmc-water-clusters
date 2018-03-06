#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
starttime=time.time()

if len(sys.argv)<4:
    print 'Usage: ./HOHOH-LeftOrRight-SharedProton.py N_size nReps descendantSteps nRepsDW'
    end

N_size=int(sys.argv[1])
nReps=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])

equilibrationSteps=500

au2wn=219474.63
nBins=51
AvePsi2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]
GatherExpectationRn=[]
Wfn=dmc.wavefunction('HOHOH', N_size)
Wfn.setNodalSurface('SharedProton',side='Both')
Destination='ResultsH3O2/'    
plotFileName=Destination+'Vref-Pop-histogram-ExcState'+str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)+'.png'
#Equilibration
initialx=Wfn.x*1.1
#initialx[:N_size/2]=Wfn.x[:N_size/2]*1.1
#initialx[N_size/2:]=Wfn.x[N_size/2:]*-1.1

for iwfn in range(nReps):
    print '   REPETITION NUMBER: ', iwfn
    print 'initial V', Wfn.molecule.V(np.array([initialx[0]]))*au2wn
    print 'initial Rn', np.average(Wfn.molecule.calcRn(np.array([initialx[0]])))
    print 'equilibrating for ', equilibrationSteps, 'steps (',equilibrationSteps*Wfn.dtau,' au)'
    v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
    inputx=equilibratedCoordinates
    
    print 'equilibrated Rn', np.average(Wfn.molecule.calcRn(equilibratedCoordinates))
    plt.figure(1)
    plt.subplot(311)
    plt.plot(np.arange(equilibrationSteps+1),np.array(v_ref_equilibration)*au2wn)
    plt.subplot(312)
    plt.plot(np.arange(equilibrationSteps+1),np.array(pop_list_equilibration))

    Rn=Wfn.molecule.calcSharedProtonDisplacement(equilibratedCoordinates)
    print 'Rn for this section was', np.average(Rn)
    descendantWeights=np.zeros((equilibratedCoordinates.shape[0]))
    for ides in range(nRepsDW):
        print 'DW Rep Number',ides,
        v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(equilibratedCoordinates,descendantSteps,initialPop=N_size, printCensus=False)
        descendantWeights=descendantWeights+descendantsTemp
    descendantWeights=descendantWeights/nRepsDW

    #Wfn.exportCoords(equilibratedCoordinates,Destination+'HOHOH-Excited-Eq-'+str(iwfn)+'.xyz',descendantWeights)

    Psi2Hist,bin_edges=np.histogram(Rn, bins=nBins, range=(-2.5,2.5),density=True,weights=descendantWeights)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
    plt.subplot(313)

    plt.plot(bin_center,Psi2Hist)
    
    print '-------------'

    Expectation_Rn=np.average(Rn,weights=descendantWeights)

    if Expectation_Rn<0:
        print '       Left',
    else:
        print '       Right',
    print Expectation_Rn
    print '       Averaged V_ref',np.average(np.array(v_ref_equilibration)[200:])*au2wn, '  cm-1',
    print 'std:',np.std(v_ref_equilibration)*au2wn, ' cm-1, unc:', (np.max(v_ref_equilibration)-np.min(v_ref_equilibration))/(2.0*np.sqrt(equilibrationSteps-200))*au2wn
    GatherExpectationRn.append(np.sum(Rn*Rn*descendantWeights)/np.sum(descendantWeights))
    averaged_vref.append(np.average(np.array(v_ref_equilibration)[200:])*au2wn)
endtime=time.time()


plt.savefig(plotFileName)
plt.clf()
print '--------   Eo   --------'
print 'averaged v_ref:',averaged_vref
print 'the average of average V_ref is',np.average(np.array(averaged_vref)), ' cm-1',
print 'standard deviation', np.std(np.array(averaged_vref)), ' cm-1'
print 'uncertainity is', (np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps))
print '--------   Rn   --------' 
print '   Average:',np.average(GatherExpectationRn),'\n   Standard Deviation:',np.std(GatherExpectationRn)
print '   Uncertainity:',(np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps))

fileR2Data=open('R2Data-ExcitedFixedNodeState.data','a')
fileR2Data.write('1   1 '+str(N_size)+'   '+str(equilibrationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
fileR2Data.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
fileR2Data.write(str(np.average(GatherExpectationRn))+'     '+str(np.std(GatherExpectationRn))+'     ')
fileR2Data.write(str((np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps)))+'\n')
fileR2Data.close()


print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'


print 'done!'
