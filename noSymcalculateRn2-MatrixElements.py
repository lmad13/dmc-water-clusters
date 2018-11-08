#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import time
starttime=time.time()
#This is for calculating <1|R^2|0>, <1|R^2|1>, <1_0|R^2|1_0> for the 3 fixed node coordinates (deltaR, Zdisp, Asym St.)


if len(sys.argv)<6:
    print 'Usage: ./thisScript.py MOLECULE={H3O2,D3O2} STATE={AsymSt,DeltaR,Ground,ZDisp} N_size nReps descendantSteps nRepsDW '
    end


stateGround='stateGround'

molecule=sys.argv[1]
state='state'+sys.argv[2]
DWstate='DW'+sys.argv[2]
N_size=int(sys.argv[3])
nReps=int(sys.argv[4])
descendantSteps=int(sys.argv[5])
nRepsDW=int(sys.argv[6])
dTau=10
Wfn=dmc.wavefunction('HOHOH', N_size)

excitedPath='data'+molecule+'/'+state+'/'+DWstate+'/'
print 'path: ', excitedPath

groundPath='data'+molecule+'/stateGround/DWGround/'
print 'Ground Path: ', groundPath

groundToExcPath='data'+molecule+'/'+stateGround+'/'+DWstate+'/'
print 'Ground to Excited State path',groundToExcPath

Wfn=dmc.wavefunction('HOHOH', N_size)
if 'D' in molecule:
    Wfn.setIsotope('fullyDeuterated')

if sys.argv[2]=='DeltaR':
    func=Wfn.molecule.calcSharedProtonDisplacement
    name='DeltaR'
elif sys.argv[2]=='AsymSt':
    func=Wfn.molecule.calcStretchAnti
    name='AntiSymmetricStretch'
elif sys.argv[2]=='ZDisp':
    func=Wfn.molecule.calcZdisplacement
    name='ZDisplacement'
else:
    print '\nINVALID CHOICE of STATE!\n'
    print 'Choose from: {AsymSt,DeltaR,Ground,ZDisp}'
    die

#for (func,name) in [(Wfn.molecule.calcSharedProtonDisplacement,'DeltaR'),
#                    (Wfn.molecule.calcStretchAnti,'AntiSymmetricStretch'),
#                    (Wfn.molecule.calcZdisplacement,'ZDisplacement')]:
                    #(Wfn.molecule.calcDipole,'Dipole')]:
print '----------',name,'---------------'
Expectation_0_R_0=[]         ##   < 0|Rn|0 >
Expectation_0_R_1=[]         ##   < 0|Rn|1 >
Expectation_0_R_1_ver2=[]         ##   < 0|Rn|1 >
Expectation_1_R_1=[]         ##   < 1|Rn|1 > 
Expectation_1_R_1=[]         ##   < 1|Rn|1 > 
Expectation_1_R_1_PRIME=[]   ##   < 1 <0/0> |Rn| 1 >
Expectation_1_R_1_PRIME_ver2=[]   ##   < 1 <0/0> |Rn| 1 >


for iwfn in range(nReps):
    #load in the three wavefunctions
    groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    groundStateCoords,dw00=Wfn.loadCoords(groundPath+groundStateWfnName)
    
    symGroundStateCoords,DW0_0=Wfn.molecule.symmetrizeCoordinates(groundStateCoords,dw00,typeOfSymmetry='none')
    
    DW0_0=DW0_0/np.sum(DW0_0)
    
    groundToExcStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-'+stateGround+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    
    groundToExcStateCoords,dw01=Wfn.loadCoords(groundToExcPath+groundToExcStateWfnName)
    
    symGgroundToExcStateCoords,DW0_1=Wfn.molecule.symmetrizeCoordinates(groundToExcStateCoords,dw01,typeOfSymmetry='none')

    #DW0_1=np.concatenate((dw01,dw01))

    DW0_1=DW0_1/np.sum(DW0_1)

    excStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    excStateCoords,dw11=Wfn.loadCoords(excitedPath+excStateWfnName)
    symExcStateCoords,DW1_1=Wfn.molecule.symmetrizeCoordinates(excStateCoords,dw11,typeOfSymmetry='regular')

    DW1_1=DW1_1/np.sum(DW1_1)

    print 'loading these',groundPath+groundStateWfnName, groundToExcPath+groundToExcStateWfnName, excitedPath+excStateWfnName

    R_0=func(symGroundStateCoords)**2
    R_1=func(symExcStateCoords)**2

    Expectation_0_R_0.append(np.sum(R_0*DW0_0)/np.sum(DW0_0))
    Expectation_1_R_1.append(np.sum(R_1*DW1_1)/np.sum(DW1_1))
    nonZeros_0=np.logical_not((DW0_0==0.0))
    sign_groundState=np.sign(func(symGroundStateCoords))
    Expectation_0_R_1.append(np.sum(R_0[nonZeros_0]*sign_groundState[nonZeros_0]*DW0_1[nonZeros_0]) / np.sqrt(np.sum(DW0_0[nonZeros_0])*np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0])))

    
    Expectation_1_R_1_PRIME.append(np.sum(R_0[nonZeros_0]*DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0])/np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0]))
    
    zeros_0=np.logical_not((DW0_0>0.0))
    print 'normalization constants, <00>:',np.sum(DW0_0),'<11>:',np.sum(DW1_1),
    print '<11PRIME>',np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/(DW0_0[nonZeros_0])/np.sum(DW0_0[nonZeros_0])),
    print '<01>:',np.sqrt(np.sum(DW0_0[nonZeros_0])*np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/(DW0_0[nonZeros_0]/np.sum(DW0_0[nonZeros_0]))))
                                                                                                                                                                                                     
    print 'the old min of DW0_0:', np.min(DW0_0[nonZeros_0]),'vs.',
    
    DW0_0_ver2=1.0*DW0_0
    DW0_0_ver2[zeros_0]=np.min(DW0_0[nonZeros_0])/2.0#np.random.uniform(low=0.000000001,high=np.min(DW0_0[nonZeros_0]))
    print np.sum(DW0_1*DW0_1/DW0_0_ver2)
    
    Expectation_0_R_1_ver2.append(np.sum(R_0*DW0_1)/np.sqrt(np.sum(DW0_0_ver2)*np.sum(DW0_1*DW0_1/DW0_0_ver2)))
    Expectation_1_R_1_PRIME_ver2.append(np.sum(R_0*DW0_1*DW0_1/DW0_0_ver2)/np.sum(DW0_1*DW0_1/DW0_0_ver2))
    
    
print 'Tau DW: ', descendantSteps, "<0|",name,"|0>          =",np.average(Expectation_0_R_0),'+/-',(np.max(Expectation_0_R_0)-np.min(Expectation_0_R_0))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<0|",name,"|1>          =",np.average(Expectation_0_R_1),'+/-',(np.max(Expectation_0_R_1)-np.min(Expectation_0_R_1))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<0|",name,"|1>(ver2)    =",np.average(Expectation_0_R_1_ver2),'+/-',(np.max(Expectation_0_R_1_ver2)-np.min(Expectation_0_R_1_ver2))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1>          =",np.average(Expectation_1_R_1),'+/-',(np.max(Expectation_1_R_1)-np.min(Expectation_1_R_1))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1> Prime    =",np.average(Expectation_1_R_1_PRIME),'+/-',(np.max(Expectation_1_R_1_PRIME)-np.min(Expectation_1_R_1_PRIME))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1> Prime v2 =",np.average(Expectation_1_R_1_PRIME_ver2),'+/-',(np.max(Expectation_1_R_1_PRIME_ver2)-np.min(Expectation_1_R_1_PRIME_ver2))/(2.0*np.sqrt(nReps))
                                                                               
fileout=open(excitedPath+'MatrixElements-squared-'+name+'.data','a')
fileout.write( molecule+'-'+state+'-'+DWstate+'   '+str(dTau)+'   '+str(N_size)+'      '+str(nReps)+'      '+str(descendantSteps)+'      '+str(nRepsDW)+'      '+
               str(np.average(Expectation_0_R_0))+' +/- '+str((np.max(Expectation_0_R_0)-np.min(Expectation_0_R_0))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Expectation_0_R_1))+' +/- '+str((np.max(Expectation_0_R_1)-np.min(Expectation_0_R_1))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Expectation_0_R_1_ver2))+' +/- '+str((np.max(Expectation_0_R_1_ver2)-np.min(Expectation_0_R_1_ver2))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Expectation_1_R_1))+' +/- '+str((np.max(Expectation_1_R_1)-np.min(Expectation_1_R_1))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Expectation_1_R_1_PRIME))+' +/- '+str((np.max(Expectation_1_R_1_PRIME)-np.min(Expectation_1_R_1_PRIME))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Expectation_1_R_1_PRIME_ver2))+' +/- '+str((np.max(Expectation_1_R_1_PRIME_ver2)-np.min(Expectation_1_R_1_PRIME_ver2))/(2.0*np.sqrt(nReps)))+'\n')
fileout.close()
print 'that took', time.time()-starttime ,'s'
