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

equilibrationSteps=500
propagationSteps=500

starttime=time.time()
au2wn=219474.63
nBins=51
AvePsi2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size)
Destination='ResultsH3O2/'
GatherExpectationRn=[]
GatherExpectationRn2=[]
GatherExpectationMagMu=[]
GatherExpectationMagMu2=[]
#Equilibration
initialx=Wfn.x*1.1

print 'initial V', Wfn.molecule.V([initialx[0]])*au2wn
print 'equilibrating for ', equilibrationSteps, 'steps (',equilibrationSteps*Wfn.dtau,' au)'
v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
inputx=equilibratedCoordinates
plt.figure(1)
plt.subplot(311)
plt.plot(np.arange(equilibrationSteps+1),np.array(v_ref_equilibration)*au2wn)
plt.subplot(312)
plt.plot(np.arange(equilibrationSteps+1),np.array(pop_list_equilibration))
#Sampling of Psi

for iwfn in range(nReps):
    print '   REPETITION NUMBER: ', iwfn
    v_ref_list,pop_list,finalCoords,d=Wfn.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)

    averaged_vref.append(np.average(np.array(v_ref_list)*au2wn))
    list_of_pop_list.append(pop_list)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps)+equilibrationSteps,np.array(v_ref_list)*au2wn)

    plt.figure(1)
    plt.subplot(312)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps)+equilibrationSteps,np.array(pop_list))

    Rn=Wfn.molecule.calcSharedProtonDisplacement(finalCoords)
    Dipole=Wfn.molecule.calcDipole(finalCoords)

    descendantWeights=np.zeros((finalCoords.shape[0]))

    for ides in range(nRepsDW):
        print 'DW Rep Number',ides,
        v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(finalCoords,descendantSteps,initialPop=N_size,printCensus=False)
        descendantWeights=descendantWeights+descendantsTemp

    descendantWeights=descendantWeights/nRepsDW

    #Wfn.exportCoords(inputx,Destination+'HOHOH-Ground-Eq-'+str(iwfn)+'.xyz',descendantWeights)
    Psi2Hist,bin_edges=np.histogram(Rn, bins=nBins, range=(-2.5,2.5),density=True,weights=descendantWeights)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
    plt.subplot(313)

    plt.plot(bin_center,Psi2Hist)
    AvePsi2Hist=AvePsi2Hist+Psi2Hist
    inputx=finalCoords
    print 'Dipole sample',Dipole[0:5]
    print 'norm',np.linalg.norm(Dipole[0:5],axis=0)
    print 'DW', descendantWeights[0:5]
    GatherExpectationRn2.append(np.sum(Rn*Rn*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationRn.append(np.sum(Rn*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationMagMu.append(np.sum(np.linalg.norm(Dipole,axis=1)*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationMagMu2.append(np.sum(np.linalg.norm(Dipole,axis=1)**2*descendantWeights)/np.sum(descendantWeights))
    print 'Rn^2:',GatherExpectationRn2[-1]
    print 'Rn:  ',GatherExpectationRn[-1]
    print 'Mu^2:',GatherExpectationMagMu2[-1]
    print 'Mu:  ',GatherExpectationMagMu[-1]

endtime=time.time()
plt.savefig(Destination+'Vref-Pop-histogram-GroundState.png')
plt.clf()

print 'averaged v_ref:',averaged_vref
print 'the average of average V_ref is',np.average(np.array(averaged_vref)), ' cm-1',
print 'standard deviation', np.std(np.array(averaged_vref)), ' cm-1'
print 'uncertainity is', (np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps))
print '--------   Rn   --------' 
print '   Average:',np.average(GatherExpectationRn),'\n   Standard Deviation:',np.std(GatherExpectationRn)
print '   Uncertainity:',(np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps))
print '--------   Rn^2  --------' 
print '   Average:',np.average(GatherExpectationRn2),'\n   Standard Deviation:',np.std(GatherExpectationRn2)
print '   Uncertainity:',(np.max(GatherExpectationRn2)-np.min(GatherExpectationRn2))/(2.0*np.sqrt(nReps))
print '--------   |Mu|  --------' 
print '   Average:',np.average(GatherExpectationMagMu),'\n   Standard Deviation:',np.std(GatherExpectationMagMu)
print '   Uncertainity:',(np.max(GatherExpectationMagMu)-np.min(GatherExpectationMagMu))/(2.0*np.sqrt(nReps))
print '--------   |Mu^2|  --------' 
print '   Average:',np.average(GatherExpectationMagMu2),'\n   Standard Deviation:',np.std(GatherExpectationMagMu2)
print '   Uncertainity:',(np.max(GatherExpectationMagMu2)-np.min(GatherExpectationMagMu2))/(2.0*np.sqrt(nReps))
#for pop_list in list_of_pop_list:
#    plt.plot(pop_list)
#plt.show()

fileOutData=open('R2Data-GroundState.data','a')
fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
fileOutData.write(str(np.average(GatherExpectationRn2))+'     '+str(np.std(GatherExpectationRn2))+'     ')
fileOutData.write(str((np.max(GatherExpectationRn2)-np.min(GatherExpectationRn2))/(2.0*np.sqrt(nReps)))+'\n')
fileOutData.close()

fileOutData=open('RData-GroundState.data','a')
fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
fileOutData.write(str(np.average(GatherExpectationRn))+'     '+str(np.std(GatherExpectationRn))+'     ')
fileOutData.write(str((np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps)))+'\n')
fileOutData.close()

fileOutData=open('MagMuData-GroundState.data','a')
fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
fileOutData.write(str(np.average(GatherExpectationMagMu))+'     '+str(np.std(GatherExpectationMagMu))+'     ')
fileOutData.write(str((np.max(GatherExpectationMagMu)-np.min(GatherExpectationMagMu))/(2.0*np.sqrt(nReps)))+'\n')
fileOutData.close()

fileOutData=open('MagMu2Data-GroundState.data','a')
fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
fileOutData.write(str(np.average(GatherExpectationMagMu2))+'     '+str(np.std(GatherExpectationMagMu2))+'     ')
fileOutData.write(str((np.max(GatherExpectationMagMu2)-np.min(GatherExpectationMagMu2))/(2.0*np.sqrt(nReps)))+'\n')
fileOutData.close()

print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'


print 'done!'
