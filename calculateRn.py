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

##Calculate <0|O|1>
nReps,nRepsDW=10,10
for tauDW in [descendantSteps,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        groundStateCoords,DW1_0=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Load-Ground-PropagatedOn-FixedNode-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn01=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords)
        Dipole01=Wfn.molecule.calcDipole(groundStateCoords)
        #print 'sample dipoles', Dipole01[0:3]
        magRn=np.absolute(Rn01)              

        magDip=np.linalg.norm(Dipole01,axis=1)
        #print 'sample magDip', magDip[0:3]
        ExpectationRn.append(np.sum(magRn*DW1_0)/np.sum(DW1_0))
        ExpectationDipole.append(np.sum(magDip*DW1_0)/np.sum(DW1_0))
        
    print 'Tau DW: ', tauDW, "<0| |Rn| |1> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Mu| |1> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))
                                                                                   

##Calculate <1_0|O|1_0>
nReps,nRepsDW=10,10
for tauDW in [descendantSteps,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    MatrixRn=[]
    MatrixDipole=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        groundStateCoords,DW00=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Ground-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        groundStateCoords_0,DW1_0=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Load-Ground-PropagatedOn-FixedNode-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords)
        Dipole0=Wfn.molecule.calcDipole(groundStateCoords)
        #Rn0_0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords_0)
        print 'same?', np.allclose(groundStateCoords, groundStateCoords_0)
        magRn=np.absolute(Rn0)                  
        magDip=np.linalg.norm(Dipole0,axis=1)
        nonZeros=np.logical_not((DW00==0.0))
        MatrixRn.append(np.sum(magRn*DW1_0)/np.sum(DW1_0))
        MatrixDipole.append(np.sum(magDip*DW1_0)/np.sum(DW1_0))
        ExpectationRn.append(np.sum(magRn[nonZeros]*(DW1_0[nonZeros]**2)/DW00[nonZeros])/np.sum(DW1_0[nonZeros]**2/DW00[nonZeros]))
        ExpectationDipole.append(np.sum(magDip[nonZeros]*DW1_0[nonZeros]**2/DW00[nonZeros])/np.sum(DW1_0[nonZeros]**2/DW00[nonZeros]))
    print 'Tau DW: ', tauDW, "<0| |Rn| |1> =",np.average(MatrixRn),'+/-',(np.max(MatrixRn)-np.min(MatrixRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Mu| |1> =",np.average(MatrixDipole),'+/-',(np.max(MatrixDipole)-np.min(MatrixDipole))/(2.0*np.sqrt(nReps))

    print 'Tau DW: ', tauDW, "<1_0| |Rn| |1_0> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1_0| |Mu| |1_0> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))

##Calculate <1|O|1>
nReps,nRepsDW=5,5

for tauDW in [descendantSteps,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        excStateCoords,DW11=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn1=Wfn.molecule.calcSharedProtonDisplacement(excStateCoords)
        Dipole1=Wfn.molecule.calcDipole(excStateCoords)
        magRn=np.absolute(Rn1)
        magDip=np.linalg.norm(Dipole1,axis=1)
        ExpectationRn.append(np.sum(magRn*DW11)/np.sum(DW11))
        ExpectationDipole.append(np.sum(magDip*DW11)/np.sum(DW11))
    print 'Tau DW: ', tauDW, "<1| |Rn| |1> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<1| |Mu| |1> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))


##Calculate <0|0|0>
nReps,nRepsDW=10,10
for tauDW in [descendantSteps,2*descendantSteps,3*descendantSteps,4*descendantSteps,5*descendantSteps]:
    ExpectationRn=[]
    ExpectationDipole=[]
    for iwfn in range(nReps):
        parameterString=str(N_size)+'-'+str(nReps)+'-'+str(tauDW)+'-'+str(nRepsDW)
        groundStateCoords,DW00=Wfn.loadCoords('Wfn-HOHOH-Tau/HOHOH-Ground-'+parameterString+'Eq-'+str(iwfn)+'.xyz')
        
        Rn0=Wfn.molecule.calcSharedProtonDisplacement(groundStateCoords)
        Dipole0=Wfn.molecule.calcDipole(groundStateCoords)
        magRn=np.absolute(Rn0)                                
        magDip=np.linalg.norm(Dipole0,axis=1)
        ExpectationRn.append(np.sum(magRn*DW00)/np.sum(DW00))
        ExpectationDipole.append(np.sum(magDip*DW00)/np.sum(DW00))
    print 'Tau DW: ', tauDW, "<0| |Rn| |0> =",np.average(ExpectationRn),'+/-',(np.max(ExpectationRn)-np.min(ExpectationRn))/(2.0*np.sqrt(nReps))
    print 'Tau DW: ', tauDW, "<0| |Mu| |0> =",np.average(ExpectationDipole),'+/-',(np.max(ExpectationDipole)-np.min(ExpectationDipole))/(2.0*np.sqrt(nReps))

    
