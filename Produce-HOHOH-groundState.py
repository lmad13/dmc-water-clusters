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


state='stateGround'
DWstate='DWGround'
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
print '   *** printing to ***    ',path+fileParameterName+'-logFile.data'
equilibrationSteps=500
propagationSteps=500

starttime=time.time()

averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size)


#Equilibration
initialx=Wfn.x*1.1

print 'equilibrating for ', equilibrationSteps, 'steps (',equilibrationSteps*Wfn.dtau,' au)'
outputFile.write('equilibrating for '+str(equilibrationSteps)+' steps ('+str(equilibrationSteps*Wfn.dtau)+' au) \n')

v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
inputx=equilibratedCoordinates
outputFile.write('From equilibration, the average energy is '+str(np.average(v_ref_equilibration[equilibrationSteps/2:])*au2wn)+' cm-1 \n')
#print use.printArray(np.array(zip(np.arange(equilibrationSteps+1),v_ref_equilibration,pop_list_equilibration)))
plt.figure(1)
plt.subplot(311)

plt.plot(np.arange(equilibrationSteps+1),np.array(v_ref_equilibration)*au2wn)
plt.subplot(312)
plt.plot(np.arange(equilibrationSteps+1),pop_list_equilibration)
np.savetxt(path+'equilibration-vref-pop-'+fileParameterName+'-array.data',np.array(zip(np.arange(equilibrationSteps+1),v_ref_equilibration,pop_list_equilibration)))
outputFile.write('save energy to '+path+'equilibration-vref-pop-'+fileParameterName+'-array.data'+'\n')

#

for iwfn in range(nReps):
    print '   REPETITION NUMBER: ', iwfn
    v_ref_list,pop_list,finalCoords,d=Wfn.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)
    np.savetxt(path+'vref-pop-'+str(iwfn)+'-'+fileParameterName+'-array.data',np.array(zip(np.arange(propagationSteps+1),v_ref_list,pop_list)))
    averaged_vref.append(np.average(np.array(v_ref_list[propagationSteps/2:])*au2wn))
    print 'average from this simulation',averaged_vref[-1]
    outputFile.write('average from simulation '+str(iwfn)+' is '+str(averaged_vref[-1])+'\n')

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
    print ''
    Wfn.exportCoords(finalCoords,path+'Wfn-'+str(iwfn)+'-'+fileParameterName+'.xyz',descendantWeights)
    Psi2Hist,bin_edges=np.histogram(Rn, bins=nBins, range=(-2.5,2.5),density=True,weights=descendantWeights)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
    plt.subplot(313)
    plt.plot(bin_center,Psi2Hist)


    Psi2R2Hist,bin_edges_r2=np.histogram(Rn**2, bins=nBins, range=(-6.5,6.5),density=True,weights=descendantWeights)
    bin_center_r2=(bin_edges_r2[:-1]+bin_edges_r2[1:])/2.0


    Psi2DipHist,bin_edges_dip=np.histogram(np.linalg.norm(Dipole,axis=1), bins=nBins, range=(0,2.5),density=True,weights=descendantWeights)
    bin_center_dip=(bin_edges_dip[:-1]+bin_edges_dip[1:])/2.0


    Psi2Dip2Hist,bin_edges_dip2=np.histogram(np.linalg.norm(Dipole,axis=1)**2, bins=nBins, range=(0,8.0),density=True,weights=descendantWeights)
    bin_center_dip2=(bin_edges_dip2[:-1]+bin_edges_dip2[1:])/2.0



    inputx=finalCoords

    

endtime=time.time()

plt.show()

print 'averaged v_ref:',averaged_vref
print 'the average of average V_ref is',np.average(np.array(averaged_vref)), ' cm-1',
print 'standard deviation', np.std(np.array(averaged_vref)), ' cm-1'
print 'uncertainity is', (np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps))
outputFile.write('averaged v_ref:'+str(averaged_vref)+'\n')
outputFile.write('the average of average V_ref is'+str(np.average(np.array(averaged_vref)))+ ' cm-1'+'\n')
outputFile.write('standard deviation'+ str(np.std(np.array(averaged_vref)))+ ' cm-1'+'\n')
outputFile.write('uncertainity is'+ str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'\n')

outputFile.close()
#####print '--------   Rn   --------' 
#####print '   Average:',np.average(GatherExpectationRn),'\n   Standard Deviation:',np.std(GatherExpectationRn)
#####print '   Uncertainity:',(np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps))
#####print '--------   Rn^2  --------' 
#####print '   Average:',np.average(GatherExpectationRn2),'\n   Standard Deviation:',np.std(GatherExpectationRn2)
#####print '   Uncertainity:',(np.max(GatherExpectationRn2)-np.min(GatherExpectationRn2))/(2.0*np.sqrt(nReps))
#####print '--------   |Mu|  --------' 
#####print '   Average:',np.average(GatherExpectationMagMu),'\n   Standard Deviation:',np.std(GatherExpectationMagMu)
#####print '   Uncertainity:',(np.max(GatherExpectationMagMu)-np.min(GatherExpectationMagMu))/(2.0*np.sqrt(nReps))
#####print '--------   |Mu^2|  --------' 
#####print '   Average:',np.average(GatherExpectationMagMu2),'\n   Standard Deviation:',np.std(GatherExpectationMagMu2)
#####print '   Uncertainity:',(np.max(GatherExpectationMagMu2)-np.min(GatherExpectationMagMu2))/(2.0*np.sqrt(nReps))
#for pop_list in list_of_pop_list:
#    plt.plot(pop_list)
#plt.show()

#####fileOutData=open('R2Data-GroundState.data','a')
#####fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
#####fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
#####fileOutData.write(str(np.average(GatherExpectationRn2))+'     '+str(np.std(GatherExpectationRn2))+'     ')
#####fileOutData.write(str((np.max(GatherExpectationRn2)-np.min(GatherExpectationRn2))/(2.0*np.sqrt(nReps)))+'\n')
#####fileOutData.close()
#####
#####fileOutData=open('RData-GroundState.data','a')
#####fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
#####fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
#####fileOutData.write(str(np.average(GatherExpectationRn))+'     '+str(np.std(GatherExpectationRn))+'     ')
#####fileOutData.write(str((np.max(GatherExpectationRn)-np.min(GatherExpectationRn))/(2.0*np.sqrt(nReps)))+'\n')
#####fileOutData.close()
#####
#####fileOutData=open('MagMuData-GroundState.data','a')
#####fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
#####fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
#####fileOutData.write(str(np.average(GatherExpectationMagMu))+'     '+str(np.std(GatherExpectationMagMu))+'     ')
#####fileOutData.write(str((np.max(GatherExpectationMagMu)-np.min(GatherExpectationMagMu))/(2.0*np.sqrt(nReps)))+'\n')
#####fileOutData.close()
#####
#####fileOutData=open('MagMu2Data-GroundState.data','a')
#####fileOutData.write('0   0 '+str(N_size)+'   '+str(propagationSteps)+'   '+str(nReps)+'   '+str(descendantSteps)+'   '+str(nRepsDW)+'      ')
#####fileOutData.write(str(np.average(averaged_vref))+'   '+str(np.std(averaged_vref))+'   '+str((np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps)))+'   ')
#####fileOutData.write(str(np.average(GatherExpectationMagMu2))+'     '+str(np.std(GatherExpectationMagMu2))+'     ')
#####fileOutData.write(str((np.max(GatherExpectationMagMu2)-np.min(GatherExpectationMagMu2))/(2.0*np.sqrt(nReps)))+'\n')
#####fileOutData.close()
#####
outputFile.write('that took '+str(endtime-starttime)+' seconds and '+str((endtime-starttime)/60.0)+' minutes \n')


print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'


print 'done!'
