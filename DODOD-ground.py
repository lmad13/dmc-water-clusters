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
Wfn.setIsotope('fullyDeuterated')
inputx=Wfn.x*1.05

Destination='ResultsD3O2/'
parameterString=str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)
plotFileName=Destination+'AsyStretch-Ground-'+parameterString

##I need to measure the ZPE first!!!
for iwfn in range(nReps):
    print '\n   REPETITION NUMBER: ', iwfn
    
    v_ref_list,pop_list,finalCoords,d=Wfn.propagate(inputx,propagationSteps,printCensus=True,initialPop=N_size)
    averaged_vref.append( np.average(v_ref_list[propagationSteps/2:]))

    plt.subplot(311)

    plt.plot(np.arange(nReps*propagationSteps,propagationSteps+nReps*propagationSteps+1),np.array(v_ref_list)*au2wn)
    plt.subplot(312)
    plt.plot(np.arange(nReps*propagationSteps,propagationSteps+nReps*propagationSteps+1),np.array(pop_list))

    for tauDW in [descendantSteps]:#,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:                    
        descendantWeights=np.zeros((finalCoords.shape[0]))
        for ides in range(nRepsDW):
            print 'DW Rep Number',ides,
            v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn.propagate(finalCoords,descendantSteps,initialPop=N_size,printCensus=False)
            descendantWeights=descendantWeights+descendantsTemp
        descendantWeights=descendantWeights/nRepsDW
        #inputx=finalCoords                        

        ASymCoord,swapped=Wfn.molecule.calcStretchAnti(finalCoords)
        HistTemp,bin_edges=np.histogram(ASymCoord,bins=nBins,range=(-1.5,1.5),density=True,
                                        weights=descendantWeights)
        bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
        plt.subplot(313)
        plt.plot(bin_center,HistTemp)#,color=linecolor[iwfn])
        InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(finalCoords,descendantWeights)
        print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
        for iname,val in zip(InternalNames,Values):
            print iname, val
        Wfn.exportCoords(finalCoords,'Wfn-DODOD-Tau/DODOD-Ground-'+parameterString+'Eq-'+str(iwfn)+'.xyz',descendantWeights)





endtime=time.time()
print 'average Energy:', np.average(averaged_vref)*au2wn, np.std(averaged_vref)*au2wn,'1/cm'

print np.array(averaged_vref)*au2wn
print 'uncertainity',((np.max(averaged_vref)-np.min(averaged_vref))/(np.sqrt(2.0)*nReps)*au2wn)
print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'

plt.savefig(plotFileName+'.png')
plt.show()
plt.clf()
