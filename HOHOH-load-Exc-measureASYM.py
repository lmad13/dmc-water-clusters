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
au2ang=0.529177249
nBins=51
AvePsi2Hist=np.zeros((nBins))
AvePsi2DipHist=np.zeros((nBins))
AvePsi2R2Hist=np.zeros((nBins))
AvePsi2Dip2Hist=np.zeros((nBins))
averaged_vref=[]
list_of_pop_list=[]

Wfn=dmc.wavefunction('HOHOH', N_size)
Destination='ResultsH3O2/'
GatherExpectationRn=[]
GatherExpectationRn2=[]
GatherExpectationMagMu=[]
GatherExpectationMagMu2=[]
GatherVarianceRn=[]
GatherExpectationOH=[]

parameterString=str(N_size)+'-'+str(nReps)+'-'+str(descendantSteps)+'-'+str(nRepsDW)
for iwfn in range(nReps):
    finalCoords,descendantWeights=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-ExcitedStretchAnti-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
    
    #print '<V_ref>', np.average(Wfn.molecule.V(finalCoords)*au2wn)
    Rn=Wfn.molecule.calcSharedProtonDisplacement(finalCoords)
    Dipole=Wfn.molecule.calcDipole(finalCoords)


    Psi2Hist,bin_edges=np.histogram(Rn, bins=nBins, range=(-5.0,5.0),density=True,weights=descendantWeights)
    bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
    plt.subplot(414)
    plt.plot(bin_center,Psi2Hist)
    AvePsi2Hist=AvePsi2Hist+Psi2Hist

    Psi2R2Hist,bin_edges_r2=np.histogram(Rn**2, bins=nBins, range=(-6.5,6.5),density=True,weights=descendantWeights)
    bin_center_r2=(bin_edges_r2[:-1]+bin_edges_r2[1:])/2.0
    AvePsi2R2Hist=AvePsi2R2Hist+Psi2R2Hist

    Psi2DipHist,bin_edges_dip=np.histogram(np.linalg.norm(Dipole,axis=1), bins=nBins, range=(0,2.5),density=True,weights=descendantWeights)
    bin_center_dip=(bin_edges_dip[:-1]+bin_edges_dip[1:])/2.0
    AvePsi2DipHist=AvePsi2DipHist+Psi2DipHist

    Psi2Dip2Hist,bin_edges_dip2=np.histogram(np.linalg.norm(Dipole,axis=1)**2, bins=nBins, range=(0,8.0),density=True,weights=descendantWeights)
    bin_center_dip2=(bin_edges_dip2[:-1]+bin_edges_dip2[1:])/2.0
    AvePsi2Dip2Hist=AvePsi2Dip2Hist+Psi2Dip2Hist

    OH1,OH2,OHSym, OHAsym=Wfn.molecule.calcOHBondLenghts(finalCoords)
    Psi2OH1,bin_edges_OH=np.histogram(OH1, bins=nBins, range=(1.5,4.5),density=True,weights=descendantWeights)
    Psi2OH2,bin_edges_OH=np.histogram(OH2, bins=nBins, range=(1.5,4.5),density=True,weights=descendantWeights)
    Psi2OHSym,bin_edges_OH_sym=np.histogram(OHSym, bins=nBins, range=(0,5.0),density=True,weights=descendantWeights)
    Psi2OHAsym,bin_edges_OH_asym=np.histogram(OHAsym, bins=nBins, range=(-2.0,2.0),density=True,weights=descendantWeights)

    bin_center_OH=(bin_edges_OH[:-1]+bin_edges_OH[1:])/2.0
    bin_center_OH_sym=(bin_edges_OH_sym[:-1]+bin_edges_OH_sym[1:])/2.0
    bin_center_OH_asym=(bin_edges_OH_asym[:-1]+bin_edges_OH_asym[1:])/2.0

    plt.subplot(411)
    plt.plot(bin_center_OH,Psi2OH1)
    plt.plot(bin_center_OH,Psi2OH2)

    plt.subplot(412)
    plt.plot(bin_center_OH_asym,Psi2OHAsym)
    plt.subplot(413)
    plt.plot(bin_center_OH_sym,Psi2OHSym)

    GatherExpectationRn2.append(np.sum(Rn*Rn*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationRn.append(np.sum(Rn*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationMagMu.append(np.sum(np.linalg.norm(Dipole,axis=1)*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationMagMu2.append(np.sum(np.linalg.norm(Dipole,axis=1)**2*descendantWeights)/np.sum(descendantWeights))
    GatherVarianceRn.append(np.sum((Rn-GatherExpectationRn[-1])**2*descendantWeights)/np.sum(descendantWeights))
    GatherExpectationOH.append([np.sum(OH1*descendantWeights)/np.sum(descendantWeights),np.sum(OH2*descendantWeights)/np.sum(descendantWeights),np.sum(OHSym*descendantWeights)/np.sum(descendantWeights),np.sum(OHAsym*descendantWeights)/np.sum(descendantWeights)])

    print 'Rn^2:',GatherExpectationRn2[-1],
    print 'Rn:  ',GatherExpectationRn[-1],
    print 'Variance in R', GatherVarianceRn[-1], GatherExpectationRn2[-1]-GatherExpectationRn[-1]**2
    print 'Mu^2:',GatherExpectationMagMu2[-1],
    print 'Mu:  ',GatherExpectationMagMu[-1]
    print 'OHs:', np.array(GatherExpectationOH[-1])*au2ang
    

endtime=time.time()

plt.show()



print 'done!'
