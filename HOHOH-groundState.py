import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
starttime=time.time()
au2wn=219474.63
nBins=51
AvePsi2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]
N_size=20000
Wfn=dmc.wavefunction('HOHOH', N_size)
Destination='ResultsH3O2/'
#Equilibration
initialx=Wfn.x*1.1
equilibrationSteps=1000
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
propagationSteps=500
descendantSteps=25
nReps=5
nRepsDW=5
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

    descendantWeights=np.zeros((finalCoords.shape[0]))

    for ides in range(nRepsDW):
        print 'DW Rep Number',ides,
        v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(finalCoords,descendantSteps,initialPop=N_size)
        descendantWeights=descendantWeights+descendantsTemp
    descendantWeights=descendantWeights/nRepsDW
    Wfn.exportCoords(inputx,Destination+'HOHOH-Ground-Eq-'+str(iwfn)+'.xyz',descendantWeights)
    Psi2Hist,bin_edges=np.histogram(Rn, bins=nBins, range=(-2.5,2.5),density=True,weights=descendantWeights)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
    plt.subplot(313)

    plt.plot(bin_center,Psi2Hist)
    AvePsi2Hist=AvePsi2Hist+Psi2Hist
    inputx=finalCoords
    
    endtime=time.time()
plt.savefig(Destination+'Vref-Pop-histogram-GroundState.png')
plt.show()

print 'averaged v_ref:',averaged_vref
print 'the average of average V_ref for the last 100 steps is',np.average(np.array(averaged_vref)), ' cm-1',
print 'standard deviation', np.std(np.array(averaged_vref)), ' cm-1'
print 'uncertainity is', (np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps))
#for pop_list in list_of_pop_list:
#    plt.plot(pop_list)
#plt.show()



print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'


print 'done!'
