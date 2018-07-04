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

func=Wfn.molecule.calcDipole
name='Dipole'


if sys.argv[2]=='DeltaR':
    signFunc=Wfn.molecule.calcSharedProtonDisplacement
    signName='DeltaR'
elif sys.argv[2]=='AsymSt':
    signFunc=Wfn.molecule.calcStretchAnti
    signName='AntiSymmetricStretch'
elif sys.argv[2]=='ZDisp':
    signFunc=Wfn.molecule.calcZdisplacement
    signName='ZDisplacement'
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
normalization_00=[]
normalization_01=[]
normalization_11=[]
normalization_11_PRIME=[]
unNormedExp_0_R_0=[]
unNormedExp_0_R_1=[]
unNormedExp_1_R_1=[]
unNormedExp_1_R_1_PRIME=[]
summationDW01=[]
for iwfn in range(nReps):
    #load in the three wavefunctions
    groundStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    groundStateCoords,dw00=Wfn.loadCoords(groundPath+groundStateWfnName)
    
    symGroundStateCoords,DW0_0=Wfn.molecule.symmetrizeCoordinates(groundStateCoords,dw00,typeOfSymmetry='regular')

    DW0_0=DW0_0/np.sum(DW0_0)

    groundToExcStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-'+stateGround+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
    
    groundToExcStateCoords,dw01=Wfn.loadCoords(groundToExcPath+groundToExcStateWfnName)

    DW0_1=np.concatenate((dw01,dw01))
    #DW0_1=DW0_1/np.sum(DW0_1)

    excStateWfnName='Wfn-'+str(iwfn)+'-'+molecule+'-'+state+'-'+DWstate+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'

    excStateCoords,dw11=Wfn.loadCoords(excitedPath+excStateWfnName)
    symExcStateCoords,DW1_1=Wfn.molecule.symmetrizeCoordinates(excStateCoords,dw11,typeOfSymmetry='regular')
    #    excStateCoords,DW1_1=Wfn.loadCoords(excitedPath+excStateWfnName)

    DW1_1=DW1_1/np.sum(DW1_1)

    
    print 'loading these',groundPath+groundStateWfnName, groundToExcPath+groundToExcStateWfnName, excitedPath+excStateWfnName

    R_0=func(symGroundStateCoords)
    R_1=func(symExcStateCoords)
    sign_groundState=np.sign(signFunc(symGroundStateCoords))
    sign_excState=np.sign(signFunc(symExcStateCoords))  #I don't actually need this... it would be squared anyway...=1

    nonZeros_0=np.logical_not((DW0_0==0.0))

    unNormedExp_0_R_1.append(np.sum(R_0[nonZeros_0]* (sign_groundState[nonZeros_0]*DW0_1[nonZeros_0])[:,None],axis=0))
    normalization_01.append( np.sqrt(np.sum(DW0_0[nonZeros_0]) * np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0])))
    normalization_00.append( np.sum(DW0_0))
    normalization_11.append(np.sum(DW1_1))
    normalization_11_PRIME.append(np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0]))
    summationDW01.append([np.sum(DW0_1*DW0_1),np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0])])
    print 'we have removed ',(np.sum(DW0_1*DW0_1)-np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]))/(np.sum(DW0_1*DW0_1))*100,'% of descendant weights on excited state surface'
    unNormedExp_0_R_0.append(np.sum(R_0*DW0_0[:,None],axis=0))
    unNormedExp_1_R_1.append(np.sum(R_1*DW1_1[:,None],axis=0))
    unNormedExp_1_R_1_PRIME.append(np.sum(R_0[nonZeros_0]*(DW0_1[nonZeros_0]**2/DW0_0[nonZeros_0])[:,None],axis=0))
    
    Expectation_0_R_0.append(np.sum(R_0*DW0_0[:,None],axis=0)/np.sum(DW0_0))
    Expectation_1_R_1.append(np.sum(R_1*(DW1_1)[:,None],axis=0)/np.sum(DW1_1))

    Expectation_0_R_1.append(np.sum(R_0[nonZeros_0]*
                                    ( sign_groundState[nonZeros_0]*
                                      DW0_1[nonZeros_0]) [:,None]
                                    , axis=0)
                             /np.sqrt(  np.sum(DW0_0[nonZeros_0]) *
                                        np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0])))

    Expectation_1_R_1_PRIME.append(np.sum(R_0[nonZeros_0]*(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0])[:,None],axis=0)/np.sum(DW0_1[nonZeros_0]*DW0_1[nonZeros_0]/DW0_0[nonZeros_0]))
    


    DW0_0_ver2=1.0*DW0_0
    zeros_0=(DW0_0==0.0)
    DW0_0_ver2[zeros_0]=np.min(DW0_0[nonZeros_0])/2.0#np.random.uniform(low=0.000000001,high=np.min(DW0_0[nonZeros_0]))

    Expectation_0_R_1_ver2.append(np.sum(R_0*(sign_groundState*DW0_1)[:,None],axis=0)/np.sqrt(np.sum(DW0_0_ver2)*np.sum(DW0_1*DW0_1/DW0_0_ver2)))
    Expectation_1_R_1_PRIME_ver2.append(np.sum(R_0*(DW0_1*DW0_1/DW0_0_ver2)[:,None],axis=0)/np.sum(DW0_1*DW0_1/DW0_0_ver2))
    

print '<0|R|0>'
print 'unnormalized,    normalizationFactor,   normalized'
for iwfn,(un,n,e) in enumerate(zip(unNormedExp_0_R_0,normalization_00,Expectation_0_R_0)):
    print iwfn,':',un,'/',n,'=',e
    print iwfn,':',np.sum(un**2),'/',n,'=',np.sum(e**2)
print '<0|R|1>'
print 'unnormalized,    normalizationFactor,   normalized'
for iwfn,(un,n,e) in enumerate(zip(unNormedExp_0_R_1,normalization_01,Expectation_0_R_1)):
    print iwfn,':',un,'/',n,'=',e
    print iwfn,':',np.sum(un**2),'/',n,'=',np.sum(e**2)
print '<1|R|1>'
print 'unnormalized,    normalizationFactor,   normalized'
for iwfn,(un,n,e) in enumerate(zip(unNormedExp_1_R_1,normalization_11,Expectation_1_R_1)):
    print iwfn,':',un,'/',n,'=',e
    print iwfn,':',np.sum(un**2),'/',n,'=',np.sum(e**2)
print '<1|R|1> PRIME'
print 'unnormalized,    normalizationFactor,   normalized'
for iwfn,(un,n,e) in enumerate(zip(unNormedExp_1_R_1_PRIME,normalization_11_PRIME,Expectation_1_R_1_PRIME)):
    print iwfn,':',un,'/',n,'=',e
    print iwfn,':',np.sum(un**2),'/',n,'=',np.sum(e**2)

print '<0|R|0>', np.average(Expectation_0_R_0,axis=0)
print '<0|R|1>', np.average(Expectation_0_R_1,axis=0)
print '<1|R|1>', np.average(Expectation_1_R_1,axis=0)
print '<1|R|1>', np.average(Expectation_1_R_1_PRIME,axis=0)

Mu2_00=np.sum(np.array(Expectation_0_R_0)**2,axis=1)
Mu2_01=np.sum(np.array(Expectation_0_R_1)**2,axis=1)
Mu2_01_ver2=np.sum(np.array(Expectation_0_R_1_ver2)**2,axis=1)
Mu2_11=np.sum(np.array(Expectation_1_R_1)**2,axis=1)
Mu2_11_PRIME=np.sum(np.array(Expectation_1_R_1_PRIME)**2,axis=1)
Mu2_11_PRIME_ver2=np.sum(np.array(Expectation_1_R_1_PRIME_ver2)**2,axis=1)
print 'Mu2_00',Mu2_00

Mu_00=np.array(Expectation_0_R_0)
Mu_01=np.array(Expectation_0_R_1)
Mu_01_ver2=np.array(Expectation_0_R_1_ver2)
Mu_11=np.array(Expectation_1_R_1)
Mu_11_PRIME=np.array(Expectation_1_R_1_PRIME)
Mu_11_PRIME_ver2=np.array(Expectation_1_R_1_PRIME_ver2)


print 'hopefully average', np.average(Mu2_00)

    
print 'Tau DW: ', descendantSteps, "<0|",name,"|0>          =",np.average(Mu2_00,axis=0),'+/-',(np.max(Mu2_00)-np.min(Mu2_00))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<0|",name,"|1>          =",np.average(Mu2_01,axis=0),'+/-',(np.max(Mu2_01)-np.min(Mu2_01))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<0|",name,"|1>(ver2)    =",np.average(Mu2_01_ver2,axis=0),'+/-',(np.max(Mu2_01_ver2)-np.min(Mu2_01_ver2))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1>          =",np.average(Mu2_11,axis=0),'+/-',(np.max(Mu2_11)-np.min(Mu2_11))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1> Prime    =",np.average(Mu2_11_PRIME,axis=0),'+/-',(np.max(Mu2_11_PRIME)-np.min(Mu2_11_PRIME))/(2.0*np.sqrt(nReps))
print 'Tau DW: ', descendantSteps, "<1|",name,"|1> Prime v2 =",np.average(Mu2_11_PRIME_ver2,axis=0),'+/-',(np.max(Mu2_11_PRIME_ver2)-np.min(Mu2_11_PRIME_ver2))/(2.0*np.sqrt(nReps))
                                                                               
fileout=open(excitedPath+'MatrixElements-squared-'+name+'.data','a')
fileout.write( molecule+'-'+state+'-'+DWstate+'   '+str(dTau)+'   '+str(N_size)+'      '+str(nReps)+'      '+str(descendantSteps)+'      '+str(nRepsDW)+'      '+
               str(np.average(Mu2_00))+' +/- '+str((np.max(Mu2_00)-np.min(Mu2_00))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Mu2_01))+' +/- '+str((np.max(Mu2_01)-np.min(Mu2_01))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Mu2_01_ver2))+' +/- '+str((np.max(Mu2_01_ver2)-np.min(Mu2_01_ver2))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Mu2_11))+' +/- '+str((np.max(Mu2_11)-np.min(Mu2_11))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Mu2_11_PRIME))+' +/- '+str((np.max(Mu2_11_PRIME)-np.min(Mu2_11_PRIME))/(2.0*np.sqrt(nReps)))+'  '+
               str(np.average(Mu2_11_PRIME_ver2))+' +/- '+str((np.max(Mu2_11_PRIME_ver2)-np.min(Mu2_11_PRIME_ver2))/(2.0*np.sqrt(nReps)))+'\n')
fileout.close()

Exp00=np.average(Mu_00,axis=0)
Exp01=np.average(Mu_01,axis=0)
Exp01_ver2=np.average(Mu_01_ver2,axis=0)
Exp11=np.average(Mu_11,axis=0)
Exp11_PRIME=np.average(Mu_11_PRIME,axis=0)
Exp11_PRIME_ver2=np.average(Mu_11_PRIME_ver2,axis=0)
                                                                               
fileout=open(excitedPath+'ComponentsOfMatrixElements-'+name+'.data','a')
fileout.write( molecule+'-'+state+'-'+DWstate+'   '+str(dTau)+'   '+str(N_size)+'      '+str(nReps)+'      '+str(descendantSteps)+'      '+str(nRepsDW)+'      <0|mu|0> '+
               str(Exp00[0])+"   "+str(Exp00[1])+"   "+str(Exp00[2])+"   "+
               '<0|mu|1> '+ str(Exp01[0])+"   "+str(Exp01[1])+"   "+str(Exp01[2])+"   "+
               '<1|mu|1>  '    + str(Exp11[0])+"   "+str(Exp11[1])+"   "+str(Exp11[2])+"   "+
               '<1|mu|1>Prime '+ str(Exp11_PRIME[0])+"   "+str(Exp11_PRIME[1])+"   "+str(Exp11_PRIME[2])+"\n")

fileout.close()
print 'that took', time.time()-starttime ,'s'

