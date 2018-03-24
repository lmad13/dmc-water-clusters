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
Psi1Psi1=np.zeros((2,nBins))
averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size,dtau=10.0)
Wfn.setNodalSurface('OHStretchAnti',side='Both')
Destination='ResultsH3O2/'
GatherExpectationRn=[]
GatherExpectationRn2=[]
GatherExpectationMagMu=[]
GatherExpectationMagMu2=[]

initialx=Wfn.x*1.10000
print 'here is where i print all of the bondlengths'
print 'atom',0 ,'and',1 ,np.average(Wfn.molecule.bondlength(initialx,atom1=0, atom2=1))
print 'atom',0 ,'and',2 ,np.average(Wfn.molecule.bondlength(initialx,atom1=0, atom2=2))
print 'atom',0 ,'and',3 ,np.average(Wfn.molecule.bondlength(initialx,atom1=0, atom2=3))
print 'atom',0 ,'and',4 ,np.average(Wfn.molecule.bondlength(initialx,atom1=0, atom2=4))
print 'atom',1 ,'and',2 ,np.average(Wfn.molecule.bondlength(initialx,atom1=1, atom2=2))
print 'atom',1 ,'and',3 ,np.average(Wfn.molecule.bondlength(initialx,atom1=1, atom2=3))
print 'atom',1 ,'and',4 ,np.average(Wfn.molecule.bondlength(initialx,atom1=1, atom2=4))
print 'atom',2 ,'and',3 ,np.average(Wfn.molecule.bondlength(initialx,atom1=2, atom2=3))
print 'atom',2 ,'and',4 ,np.average(Wfn.molecule.bondlength(initialx,atom1=2, atom2=4))
print 'atom',3 ,'and',4 ,np.average(Wfn.molecule.bondlength(initialx,atom1=3, atom2=4))

print 

initialxR=Wfn.exchange(initialx[:N_size/2],[(1,3),(2,4)])*1.0000

ASx,listOfSwapped=Wfn.molecule.calcStretchAnti(initialx)
ASr,listOfSwapped=Wfn.molecule.calcStretchAnti(initialxR)
print 'first antisymcoordinate', np.average(ASx),np.average(ASr)

#initialxL=Wfn.x[:N_size]*1.000
print initialx.shape,initialxR.shape
initialx[:N_size/2]=initialxR*1.00000

#print 'initial antisym coordinate',np.average(Wfn.molecule.calcStretchAnti(initialx))

#print 'initial V', Wfn.molecule.V([initialx[0]])*au2wn
#print 'equilibrating for ', equilibrationSteps, 'steps (',equilibrationSteps*Wfn.dtau,' au)'
#v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,equilibrationSteps,printCensus=True,initialPop=N_size)
inputx=initialx

parameterString=str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)
plotFileName=Destination+'Vref-Pop-histogram-Excited-StretchAnti-'+parameterString
linecolor=['#ffe6e6', '#e6e6ff',
'#ff9999', '#9999ff',
'#ff4d4d', '#4d4dff',
'#ff0000', '#0000ff',
'#b30000', '#0000b3',
'#660000', '#000066',
'#330000', '#000033',
'#1a0000', '#00001a',]
plt.figure(1)

#Sampling of Psi

for iwfn in range(nReps):
    print '\n   REPETITION NUMBER: ', iwfn
    v_ref_list,pop_list,finalCoords,d=Wfn.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)
    Psi1Psi1=np.zeros((2,nBins))
    averaged_vref.append(np.average(np.array(v_ref_list[propagationSteps/2:])*au2wn))
    print 'E_ref[Nsteps/2:]',averaged_vref[-1],'-',6604,'=',averaged_vref[-1]-6604
    list_of_pop_list.append(pop_list)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps),np.array(v_ref_list)*au2wn)

    plt.figure(1)
    plt.subplot(312)
    plt.plot(np.arange(iwfn*propagationSteps-1,(iwfn+1)*propagationSteps),np.array(pop_list))

    #Rn=Wfn.molecule.calcSharedProtonDisplacement(finalCoords)
    #Dipole=Wfn.molecule.calcDipole(finalCoords)
    itau=0
    for tauDW in [descendantSteps]:#,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
        descendantWeights=np.zeros((finalCoords.shape[0]))
        for ides in range(nRepsDW):
            print 'DW Rep Number',ides,
            v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(finalCoords,descendantSteps,initialPop=N_size,printCensus=False)
            descendantWeights=descendantWeights+descendantsTemp
        descendantWeights=descendantWeights/nRepsDW
        ##inputx=finalCoords
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        Wfn.exportCoords(finalCoords,'Wfn-HOHOH-Tau/HOHOH-ExcitedStretchAnti-noRecrossing'+parameterString+'Eq-'+str(iwfn)+'.xyz',descendantWeights)
        
        ASymCoord,swapped=Wfn.molecule.calcStretchAnti(finalCoords)
        print ASymCoord.shape, descendantWeights.shape
        HistTemp,bin_edges=np.histogram(ASymCoord,bins=nBins,range=(-1.5,1.5),density=True,
                                          weights=descendantWeights)

        Psi1Psi1[itau]=Psi1Psi1[itau]+HistTemp
        itau=itau+1

    plt.subplot(313)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0

    #plt.plot(bin_center,Psi1Psi1[0])#,color=linecolor[iwfn])
    plt.plot(bin_center,Psi1Psi1[0])#,color=linecolor[iwfn])
endtime=time.time()



plt.savefig(plotFileName+'.png')
plt.show()
plt.clf()

print 'averaged_vref',averaged_vref,np.average(averaged_vref),np.std(averaged_vref)
print 'averaged_vref',np.array(averaged_vref)-6604,np.average(averaged_vref)-6604,np.std(averaged_vref)
print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'

print 'done!'
