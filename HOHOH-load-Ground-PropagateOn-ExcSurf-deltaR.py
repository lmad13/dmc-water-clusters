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
AvePsi2DipHist=np.zeros((nBins))
AvePsi2R2Hist=np.zeros((nBins))
AvePsi2Dip2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]

Wfn_Zdisp=dmc.wavefunction('HOHOH', N_size)
Wfn_Zdisp.setNodalSurface('Z-displacement',side='Both')

Wfn_deltaR=dmc.wavefunction('HOHOH', N_size)
Wfn_deltaR.setNodalSurface('SharedProton',side='Both')

Destination='Wfn-HOHOH-Tau/'
GatherExpectationRn=[]
GatherExpectationRn2=[]
GatherExpectationMagMu=[]
GatherExpectationMagMu2=[]



WfnGround=dmc.wavefunction('HOHOH', N_size)
initialx=WfnGround.x*1.1
#print 'initial V', Wfn.molecule.V([initialx[0]])*au2wn
print 'equilibrating  on Ground for ', equilibrationSteps, 'steps (',equilibrationSteps*WfnGround.dtau,' au)'
v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=WfnGround.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
inputx=equilibratedCoordinates*1.000000

parameterString=str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)
plotFileName=Destination+'Vref-Pop-histogram-Excited-deltaR'+parameterString
plt.figure(1)
plt.subplot(311)
plt.plot(np.arange(equilibrationSteps+1),np.array(v_ref_equilibration)*au2wn,color='black')
plt.subplot(312)
plt.plot(np.arange(equilibrationSteps+1),np.array(pop_list_equilibration),color='black')
RnCoord,swapped=Wfn_deltaR.molecule.calcRn(equilibratedCoordinates)

descendantWeights=np.zeros((equilibratedCoordinates.shape[0]))
for ides in range(nRepsDW):
    print 'DW Rep Number',ides,
    v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=WfnGround.propagate(equilibratedCoordinates,descendantSteps,initialPop=N_size,printCensus=False)
    descendantWeights=descendantWeights+descendantsTemp
descendantWeights=descendantWeights/nRepsDW

WfnGround.exportCoords(equilibratedCoordinates,'Wfn-HOHOH-Tau/HOHOH-Ground-'+parameterString+'Eq-0.xyz',descendantWeights)

HistTemp,bin_edges=np.histogram(RnCoord,bins=nBins,range=(-1.5,1.5),density=True )
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.subplot(313)
plt.plot(bin_center,HistTemp,color='black')#,color=linecolor[iwfn])



RnCoordsHist=np.zeros((nBins,nBins))

#Sampling of Psi

for iwfn in range(nReps):
    print '\n   REPETITION NUMBER: ', iwfn
    v_ref_list,pop_list,finalCoords,d=Wfn_deltaR.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)

    averaged_vref.append(np.average(np.array(v_ref_list[propagationSteps/2:])*au2wn))
    print 'average v_ref from this sim is ', averaged_vref[-1],'1/cm'
    list_of_pop_list.append(pop_list)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(np.arange(equilibrationSteps+iwfn*propagationSteps-1,equilibrationSteps+(iwfn+1)*propagationSteps),np.array(v_ref_list)*au2wn)
    
    plt.figure(1)
    plt.subplot(312)
    plt.plot(np.arange(equilibrationSteps+iwfn*propagationSteps-1,equilibrationSteps+(iwfn+1)*propagationSteps),np.array(pop_list))

    #Rn=Wfn_deltaR.molecule.calcSharedProtonDisplacement(finalCoords)
    #Dipole=Wfn_deltaR.molecule.calcDipole(finalCoords)

    for tauDW in [descendantSteps]:#,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
        descendantWeights=np.zeros((finalCoords.shape[0]))
        for ides in range(nRepsDW):
            print 'DW Rep Number',ides,
            v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn_deltaR.propagate(finalCoords,tauDW,initialPop=N_size,printCensus=False)
            descendantWeights=descendantWeights+descendantsTemp
        descendantWeights=descendantWeights/nRepsDW
        ##inputx=finalCoords
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        Wfn_deltaR.exportCoords(finalCoords,'Wfn-HOHOH-Tau/HOHOH-Excited-DeltaR-'+parameterString+'Eq-'+str(iwfn)+'.xyz',descendantWeights)
        RnCoord,swapped=Wfn_deltaR.molecule.calcRn(finalCoords)
        HistTemp,bin_edges=np.histogram(RnCoord,bins=nBins,range=(-1.5,1.5),density=True,
                                        weights=descendantWeights)

        bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
        plt.subplot(313)
        plt.plot(bin_center,HistTemp)#,color=linecolor[iwfn])
        print 'delta R info:', np.min(RnCoord),np.max(RnCoord),np.average(RnCoord),np.average(RnCoord,weights=descendantWeights),np.std(RnCoord)
        deltaRCoord,swapped2=Wfn_deltaR.molecule.calcRn(finalCoords)
        #histogram2d
        temp2Dhist,xedges,yedges=np.histogram2d(RnCoord,deltaRCoord,bins=nBins,range=[[-1.5,1.5],[-1.5,1.5]],weights=descendantWeights)

        RnCoordsHist=RnCoordsHist+temp2Dhist

        print 'Rn Min',np.min(RnCoord), 'Rn Max', np.max(RnCoord),'deltaRMin', np.min(deltaRCoord),'deltaRMax', np.max(deltaRCoord)
endtime=time.time()
plt.figure(2)
plt.imshow(RnCoordsHist,extent=[-1.5,1.5,-1.5,1.5])

print '\n   Average Energy:', np.average(averaged_vref), np.std(averaged_vref)

print np.array(averaged_vref)
print 'uncertainity',((np.max(averaged_vref)-np.min(averaged_vref))/(np.sqrt(2.0)*nReps))
print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'
plt.figure(1)
plt.savefig(plotFileName+'.png')
plt.show()
plt.clf()




plt.savefig(plotFileName+'.png')
plt.clf()




print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'

print 'done!'
