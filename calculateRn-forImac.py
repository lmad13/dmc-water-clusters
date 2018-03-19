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
Wfn=dmc.wavefunction('HOHOH', N_size)
nBins=51
Psi0Psi1=np.zeros((nBins))
Psi1Psi1=np.zeros((nBins))
Psi1Psi1from0=np.zeros((nBins))
Psi0Psi0=np.zeros((nBins))
##Calculate <0|O|1>
##Calculate <1_0|O|1_0>
nReps,nRepsDW=10,10
for tauDW in [descendantSteps]:#,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    ExpectationSymStretch=[]
    MatrixRn=[]
    MatrixDipole=[]
    MatrixSymStretch=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        groundStateCoords,DW00=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Ground-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        groundStateCoords_0,DW1_0=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Load-Ground-PropagatedOn-FixedNode-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords)
        Dipole0=Wfn.molecule.calcDipole(groundStateCoords)
        #Rn0_0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords_0)
        #print 'same?', np.allclose(groundStateCoords, groundStateCoords_0)
        magRn=np.absolute(Rn0)                  
        magDip=np.linalg.norm(Dipole0,axis=1)
        SymStretch=Wfn.molecule.calcStretchSym(groundStateCoords)
        AsymStretch=Wfn.molecule.calcStretchAntiIn(groundStateCoords)
        nonZeros=np.logical_not((DW00==0.0))

        MatrixRn.append(np.sum(magRn*DW1_0)/np.sum(DW1_0))
        MatrixDipole.append(np.sum(magDip*DW1_0)/np.sum(DW1_0))
        MatrixSymStretch.append(np.sum(SymStretch*DW1_0)/np.sum(DW1_0))

        ExpectationRn.append(np.sum(magRn[nonZeros]*(DW1_0[nonZeros]**2)/DW00[nonZeros])/np.sum(DW1_0[nonZeros]**2/DW00[nonZeros]))
        ExpectationDipole.append(np.sum(magDip[nonZeros]*DW1_0[nonZeros]**2/DW00[nonZeros])/np.sum(DW1_0[nonZeros]**2/DW00[nonZeros]))
        ExpectationSymStretch.append(np.sum(SymStretch[nonZeros]*DW1_0[nonZeros]**2/DW00[nonZeros])/np.sum(DW1_0[nonZeros]**2/DW00[nonZeros]))

        RnHisttemp,bin_edges=np.histogram(magRn[nonZeros],bins=nBins,range=(0,2.5),density=True, 
                                          weights=(DW1_0[nonZeros]**2)/DW00[nonZeros])
        Psi1Psi1from0=Psi1Psi1from0+RnHisttemp

        RnHisttemp,bin_edges=np.histogram(magRn,bins=nBins,range=(0,2.5),density=True, weights=DW1_0)
        Psi0Psi1=Psi0Psi1+RnHisttemp

        
    print 'Tau DW: ', tauDW, "<0| |Sym| |1> =",np.average(MatrixSymStretch),'+/-',(np.max(MatrixSymStretch)-np.min(MatrixSymStretch))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Rn| |1> =",np.average(MatrixRn),'+/-',(np.max(MatrixRn)-np.min(MatrixRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Mu| |1> =",np.average(MatrixDipole),'+/-',(np.max(MatrixDipole)-np.min(MatrixDipole))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1_0| |Sym| |1_0> =",np.average(ExpectationSymStretch),'+/-',(np.max(ExpectationSymStretch)-np.min(ExpectationSymStretch))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1_0| |Rn| |1_0> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1_0| |Mu| |1_0> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))

##Calculate <1|O|1>
nReps,nRepsDW=5,5

for tauDW in [descendantSteps]:#,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    ExpectationSymStretch=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        excStateCoords,DW11=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn1=Wfn.molecule.calcSharedProtonDisplacement(excStateCoords)
        Dipole1=Wfn.molecule.calcDipole(excStateCoords)
        magRn=np.absolute(Rn1)
        magDip=np.linalg.norm(Dipole1,axis=1)
        ExpectationRn.append(np.sum(magRn*DW11)/np.sum(DW11))
        ExpectationDipole.append(np.sum(magDip*DW11)/np.sum(DW11))
        AsymStretch=Wfn.molecule.calcStretchAntiIn(excStateCoords)
        SymStretch=Wfn.molecule.calcStretchSym(excStateCoords)
        ExpectationSymStretch.append(np.sum(SymStretch*DW11)/np.sum(DW11))

        RnHisttemp,bin_edges=np.histogram(magRn,bins=nBins,range=(0,2.5),density=True, weights=DW11)
        Psi1Psi1=Psi1Psi1+RnHisttemp

    print 'Tau DW: ', tauDW, "<1| |Sym| |1> =",np.average(ExpectationSymStretch),'+/-',(np.max(ExpectationSymStretch)-np.min(ExpectationSymStretch))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1| |Rn| |1> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1| |Mu| |1> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))


##Calculate <0|0|0>
nReps,nRepsDW=10,10
for tauDW in [descendantSteps]:#,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    ExpectationSymStretch=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        groundStateCoords,DW00=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Ground-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords)
        Dipole0=Wfn.molecule.calcDipole(groundStateCoords)
        magRn=np.absolute(Rn0)                                
        magDip=np.linalg.norm(Dipole0,axis=1)
        ExpectationRn.append(np.sum(magRn*DW00)/np.sum(DW00))
        ExpectationDipole.append(np.sum(magDip*DW00)/np.sum(DW00))

        AsymStretch=Wfn.molecule.calcStretchAntiIn(groundStateCoords)
        SymStretch=Wfn.molecule.calcStretchSym(groundStateCoords)
        ExpectationSymStretch.append(np.sum(SymStretch*DW00)/np.sum(DW00))

        RnHisttemp,bin_edges=np.histogram(magRn,bins=nBins,range=(0,2.5),density=True, weights=DW00)
        Psi0Psi0=Psi0Psi0+RnHisttemp
    print 'Tau DW: ', tauDW, "<0| |Sym| |0> =",np.average(ExpectationSymStretch),'+/-',(np.max(ExpectationSymStretch)-np.min(ExpectationSymStretch))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Rn| |0> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Mu| |0> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))

    
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0


Psi0Psi0=Psi0Psi0/10.0
Psi0Psi1=Psi0Psi1/10.0
Psi1Psi1from0=Psi1Psi1from0/10.0
Psi1Psi1=Psi1Psi1/5.0
plt.figure(figsize=(12,10))
plt.plot(bin_center,Psi0Psi0/np.sum(Psi0Psi0),color='red',linewidth=2,label='Ground State')
plt.plot(bin_center,Psi0Psi1/np.sum(Psi0Psi1), color='blue',linewidth=2,label='Projection of Ground State on Excited State')
plt.plot(bin_center,Psi1Psi1/np.sum(Psi1Psi1), color='black',linewidth=2,label='Excited State')
plt.plot(bin_center,Psi1Psi1from0/np.sum(Psi1Psi1from0),color='magenta',linewidth=2,label='Excited State Propagated from Ground State')

plt.plot(bin_center,np.sqrt(Psi1Psi1)*np.sqrt(Psi0Psi0)/np.sum(np.sqrt(Psi1Psi1)*np.sqrt(Psi0Psi0)),color='green', linewidth=2.0,label='(Ground State x Excited State)^(1/2)')
plt.ylabel('Psi^2')
plt.xlabel('|R_n| (bohr)')


plt.legend()
plt.savefig('differentMagRn.png')
plt.show()
