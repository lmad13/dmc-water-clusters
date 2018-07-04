#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import time
import usefulFunctions as use
import CalculateSpectrum

starttime=time.time()
#This is for calculating <1|Mu^2|0> for  Zdisp and  Asym St. for exactly one wavefunction (for now)


if len(sys.argv)<4:
    print 'Usage: ./thisScript.py  N_size WfnNumber descendantSteps nRepsDW '
    end


stateGround='stateGround'

molecule='H3O2'
GroundState='stateGround'
DWGroundState='DWGround'
DWAsymState='DWAsymSt'
DWZDispState='DWZDisp'

N_size=int(sys.argv[1])
nWfn=int(sys.argv[2])
descendantSteps=int(sys.argv[3])
nRepsDW=int(sys.argv[4])

dTau=10
Wfn=dmc.wavefunction('HOHOH', N_size)

#First: Ground state
groundPath='data'+molecule+'/'+GroundState+'/'+DWGroundState+'/'
print 'Ground Path: ', groundPath

groundStateWfnName='Wfn-'+str(nWfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
groundStateCoords,dw00=Wfn.loadCoords(groundPath+groundStateWfnName)
symGroundStateCoords,DW0_0=Wfn.molecule.symmetrizeCoordinates(groundStateCoords,dw00,typeOfSymmetry='regular')

symEckGroundStateCoords=Wfn.molecule.eckartRotate(symGroundStateCoords)
#Wfn.exportCoords(symEckGroundStateCoords,'checkEck-compareFNtoHOA.xyz',DW0_0)
nonZeros_0=np.logical_not((DW0_0==0.0))

#Second Asym St. Fixed Node
groundToExcPath='data'+molecule+'/'+GroundState+'/'+DWAsymState+'/'
print 'path: ',groundToExcPath
groundToExcStateWfnName='Wfn-'+str(nWfn)+'-'+molecule+'-'+stateGround+DWAsymState+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
groundToAsymStStateCoords,dw01=Wfn.loadCoords(groundToExcPath+groundToExcStateWfnName)

DW0_1_Asym=np.concatenate((dw01,dw01))
#DW0_1_Asym=dw01*1.0

DipoleMoment_0=Wfn.molecule.calcDipole(symEckGroundStateCoords,eckartRotated=True)

Asym_coord=Wfn.molecule.calcStretchAnti(symEckGroundStateCoords)
sign_groundState_Asym=np.sign(Asym_coord)

unNormedMuTerms_0_R_1_Asym = DipoleMoment_0[nonZeros_0]* (sign_groundState_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0])[:,None]

normalization_01_Asym = np.sqrt(np.sum(DW0_0[nonZeros_0]) * np.sum(DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]/DW0_0[nonZeros_0]))
print '    ***\n',np.sum(DW0_1_Asym**2),',',np.sum(DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]),'% diff', 100*(np.sum(DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0])-np.sum(DW0_1_Asym**2))/np.sum(DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]),'%\n***'


Amplitude_1_Asym= sign_groundState_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]
Amplitude_11_Asym= DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]/DW0_0[nonZeros_0]

#Second ZDisp Fixed Node
groundToExcPath='data'+molecule+'/'+GroundState+'/'+DWZDispState+'/'
print 'path: ',groundToExcPath
groundToExcStateWfnName='Wfn-'+str(nWfn)+'-'+molecule+'-'+stateGround+DWZDispState+'-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.xyz'
groundToZDispStStateCoords,dw01=Wfn.loadCoords(groundToExcPath+groundToExcStateWfnName)
DW0_1_ZDisp=np.concatenate((dw01,dw01))
#DW0_1_ZDisp=dw01*1.0

DipoleMoment_0=Wfn.molecule.calcDipole(symEckGroundStateCoords,eckartRotated=True)


ZDisp_coord=Wfn.molecule.calcZdisplacement(symEckGroundStateCoords)
sign_groundState_ZDisp=np.sign(ZDisp_coord)
Amplitude_1_ZDisp= sign_groundState_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]
Amplitude_11_ZDisp= DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]/DW0_0[nonZeros_0]

unNormedMuTerms_0_R_1_ZDisp = DipoleMoment_0[nonZeros_0]* (sign_groundState_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0])[:,None]
                            
normalization_01_ZDisp = np.sqrt(np.sum(DW0_0[nonZeros_0]) * np.sum(DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]/DW0_0[nonZeros_0]))

print '    ***\n',np.sum(DW0_1_ZDisp**2),',',np.sum(DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]),'% diff',100*(np.sum(DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0])-np.sum(DW0_1_ZDisp**2))/np.sum(DW0_1_ZDisp**2),'% \n***'

MuTerms_0_R_1_ZDisp=((DipoleMoment_0[nonZeros_0]*(sign_groundState_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0])[:,None]) 
               / np.sqrt(np.sum(DW0_0[nonZeros_0])*
                        np.sum(DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]/DW0_0[nonZeros_0])))




#Now determine the Asym St via Hermite polynomial expansion
HOASpectrum=CalculateSpectrum.HarmonicApproxSpectrum(Wfn,symEckGroundStateCoords,DW0_0,path=groundPath)

GfileName='TheGMatrix-symmetrized-C2H-'+str(nWfn)+'-'+molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)+'.gmat'

fundamentalEnergies,fundamentalIntensities, combinationBandEnergies,combinationBandIntensities=HOASpectrum.calculateSpectrum(symEckGroundStateCoords,DW0_0,groundPath+GfileName,silent=False)

for e,i in zip(fundamentalEnergies,fundamentalIntensities):
    print 'mode ',e,i


Q_Asym,Mu_Asym,normQAsym,Q_ZDisp,Mu_ZDisp,normQZDisp,normQ00,allQ=HOASpectrum.calculateHOAWfnAmplitude(symEckGroundStateCoords[nonZeros_0],DW0_0[nonZeros_0],groundPath+GfileName,silent=False)


normalization_01_Asym_array=np.zeros(Amplitude_1_Asym.shape)+normalization_01_Asym
normalization_01_ZDisp_array=np.zeros(Amplitude_1_Asym.shape)+normalization_01_ZDisp
normalization_00_array=np.zeros(Amplitude_1_Asym.shape)+np.sqrt(np.sum(DW0_0))
normalization_Asym_array=np.zeros(Amplitude_1_Asym.shape)+np.sqrt(np.sum(DW0_1_Asym[nonZeros_0]*DW0_1_Asym[nonZeros_0]/DW0_0[nonZeros_0]))
normalization_ZDisp_array=np.zeros(Amplitude_1_Asym.shape)+np.sqrt(np.sum(DW0_1_ZDisp[nonZeros_0]*DW0_1_ZDisp[nonZeros_0]/DW0_0[nonZeros_0]))

normQAsym_array=np.zeros(Amplitude_1_Asym.shape)+normQAsym
normQZDisp_array=np.zeros(Amplitude_1_Asym.shape)+normQZDisp
normQ00_array=np.zeros(Amplitude_1_Asym.shape)+normQ00
fileNameString=molecule+'-stateGround-DWGround-dt'+str(dTau)+'-nWalk'+str(N_size)+'-nT'+str(descendantSteps)+'-nDW'+str(nRepsDW)
np.savetxt('AmplitudeComparison-Asym'+fileNameString+'.data',
           zip(Asym_coord[nonZeros_0],Amplitude_1_Asym,normalization_01_Asym_array,normalization_Asym_array,normalization_00_array,Q_Asym,Q_Asym*DW0_0[nonZeros_0], normQAsym_array,normQ00_array))

np.savetxt('AmplitudeComparison-ZDisp'+fileNameString+'.data',
           zip(ZDisp_coord[nonZeros_0],Amplitude_1_ZDisp,normalization_01_ZDisp_array,normalization_ZDisp_array,normalization_00_array,Q_ZDisp,Q_ZDisp*DW0_0[nonZeros_0], normQZDisp_array,normQ00_array))



np.savetxt('whichQ.data',np.concatenate((allQ,ZDisp_coord[nonZeros_0][:,None]),axis=1))


print 'so which of these two definitions is actually normalized?'
print np.sum(Amplitude_1_Asym*Asym_coord[nonZeros_0])/normalization_01_Asym, 'or', np.sum(Q_Asym**2*DW0_0[nonZeros_0])/(normQAsym*normQ00)
print np.sum(Amplitude_1_ZDisp*ZDisp_coord[nonZeros_0])/normalization_01_ZDisp, 'or', np.sum(Q_ZDisp**2*DW0_0[nonZeros_0])/(normQZDisp*normQ00)

print 'ratios: <1|x|0> FN:',(np.sum(Amplitude_1_Asym*Asym_coord[nonZeros_0])/normalization_01_Asym)/(np.sum(Amplitude_1_ZDisp*ZDisp_coord[nonZeros_0])/normalization_01_ZDisp)
print 'ratios:<1|x|0> HOA:',(np.sum(Q_Asym**2*DW0_0[nonZeros_0])/(normQAsym*normQ00))/(np.sum(Q_ZDisp**2*DW0_0[nonZeros_0])/(normQZDisp*normQ00))

print np.sum(unNormedMuTerms_0_R_1_Asym,axis=0)/normalization_01_Asym, 'or', np.sum(DipoleMoment_0[nonZeros_0]*(Q_Asym*DW0_0[nonZeros_0])[:,None],axis=0)/(normQAsym*normQ00),'or',np.sum(Mu_Asym,axis=0)/(normQAsym*normQ00)

print np.sum(unNormedMuTerms_0_R_1_ZDisp,axis=0)/normalization_01_ZDisp, 'or', np.sum(DipoleMoment_0[nonZeros_0]*(Q_ZDisp*DW0_0[nonZeros_0])[:,None],axis=0)/(normQZDisp*normQ00),'or',np.sum(Mu_ZDisp,axis=0)/(normQZDisp*normQ00)



print 'magnitudes mu^2 ZDisp',  np.sum( (np.sum(unNormedMuTerms_0_R_1_ZDisp,axis=0)/normalization_01_ZDisp)**2), 'or', np.sum( (np.sum(DipoleMoment_0[nonZeros_0]*(Q_ZDisp*DW0_0[nonZeros_0])[:,None],axis=0)/(normQZDisp*normQ00))**2),'or',np.sum( (np.sum(Mu_ZDisp,axis=0)/(normQZDisp*normQ00))**2)

print 'magnitudes mu^2 Asym',  np.sum( (np.sum(unNormedMuTerms_0_R_1_Asym,axis=0)/normalization_01_Asym)**2), 'or', np.sum( (np.sum(DipoleMoment_0[nonZeros_0]*(Q_Asym*DW0_0[nonZeros_0])[:,None],axis=0)/(normQAsym*normQ00))**2),'or',np.sum( (np.sum(Mu_Asym,axis=0)/(normQAsym*normQ00))**2)

print 'glorious ratios of |mu^2|  FN:', (np.sum( (np.sum(unNormedMuTerms_0_R_1_ZDisp,axis=0)/normalization_01_ZDisp)**2))/( np.sum( (np.sum(unNormedMuTerms_0_R_1_Asym,axis=0)/normalization_01_Asym)**2))

print 'glorious ratios of |mu^2| HOA:',(np.sum( np.sum(DipoleMoment_0[nonZeros_0]*(Q_ZDisp*DW0_0[nonZeros_0])[:,None],axis=0)/(normQZDisp*normQ00)**2))/(np.sum( np.sum(DipoleMoment_0[nonZeros_0]*(Q_Asym*DW0_0[nonZeros_0])[:,None],axis=0)/(normQAsym*normQ00)**2))

print 'glorious ratios of |mu^2| HOA_straight up:',(np.sum( np.sum(Mu_ZDisp,axis=0)/(normQZDisp*normQ00)**2))/(np.sum( np.sum(Mu_Asym,axis=0)/(normQAsym*normQ00)**2))

print 'Intensities: HOA:', np.sum( np.sum(Mu_ZDisp**2,axis=0)/((normQZDisp*normQ00)**2)),np.sum( np.sum(Mu_Asym**2,axis=0)/((normQAsym*normQ00)**2))

NBins=51
histFN_ZDisp,bin_edges=np.histogram(ZDisp_coord[nonZeros_0],bins=NBins, normed=False, weights=Amplitude_1_ZDisp/normalization_01_ZDisp)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-ZDisp.hist',zip(bin_center, histFN_ZDisp))

histFN_Asym,bin_edges=np.histogram(Asym_coord[nonZeros_0],bins=NBins, normed=False, weights=Amplitude_1_Asym/normalization_01_Asym)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-Asym.hist',zip(bin_center, histFN_Asym))

histHOA_Asym,bin_edges=np.histogram(Q_Asym,bins=NBins, normed=False, weights=Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'-Asym.hist',zip(bin_center, histHOA_Asym))

histHOA_ZDisp,bin_edges=np.histogram(Q_ZDisp,bins=NBins, normed=False, weights=Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'-ZDisp.hist',zip(bin_center, histHOA_ZDisp))


histFN_MuXAmp_Asym,bin_edges=np.histogram(Asym_coord[nonZeros_0],bins=NBins, normed=False, weights=(DipoleMoment_0[:,0]*np.sign(Asym_coord)*DW0_1_Asym)[nonZeros_0]/normalization_01_Asym)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-MuXAmp-Asym.hist',zip(bin_center, histFN_MuXAmp_Asym))

histHOA_MuXAmp_Asym,bin_edges=np.histogram(Q_Asym,bins=NBins, normed=False, weights=DipoleMoment_0[nonZeros_0][:,0]*Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'-MuXAmp-Asym.hist',zip(bin_center, histHOA_MuXAmp_Asym))

histHOA_MuXAmp_Asym_ver2,bin_edges=np.histogram(Asym_coord[nonZeros_0],bins=NBins, normed=False, weights=DipoleMoment_0[nonZeros_0][:,0]*Asym_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(Asym_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'-MuXAmp-Asym_ver2.hist',zip(bin_center, histHOA_MuXAmp_Asym_ver2))


histFN_MuXAmp_ZDisp,bin_edges=np.histogram(ZDisp_coord[nonZeros_0],bins=NBins, normed=False, weights=(DipoleMoment_0[:,0]*np.sign(ZDisp_coord)*DW0_1_ZDisp/normalization_01_ZDisp)[nonZeros_0])
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
print 'normalization FN ZDisp', np.sum(histFN_MuXAmp_ZDisp**2)
np.savetxt('FN-'+fileNameString+'-MuXAmp-ZDisp.hist',zip(bin_center, histFN_MuXAmp_ZDisp))

histHOA_MuXAmp_ZDisp,bin_edges=np.histogram(Q_ZDisp,bins=NBins, normed=False, weights=DipoleMoment_0[nonZeros_0][:,0]*Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
print 'normalization HOA Q ZDisp', np.sum(histHOA_MuXAmp_ZDisp**2)
np.savetxt('HOA-'+fileNameString+'-MuXAmp-ZDisp.hist',zip(bin_center, histHOA_MuXAmp_ZDisp))

histHOA_MuXAmp_ZDisp_ver2,bin_edges=np.histogram(ZDisp_coord[nonZeros_0],bins=NBins, normed=False, weights=DipoleMoment_0[nonZeros_0][:,0]*ZDisp_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(ZDisp_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
print 'normalization HOA ver2 ZDisp', np.sum(histHOA_MuXAmp_ZDisp_ver2**2)
np.savetxt('HOA-'+fileNameString+'-MuXAmp-ZDisp_ver2.hist',zip(bin_center, histHOA_MuXAmp_ZDisp_ver2))


histFN_DipX_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=Amplitude_1_Asym/normalization_01_Asym)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleX-Asym.hist',zip(bin_center, histFN_DipX_Asym))

histHOA_DipX_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleX-Asym.hist',zip(bin_center, histHOA_DipX_Asym))

histHOA_ver2_DipX_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=Asym_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(Asym_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleX-Asym-ver2.hist',zip(bin_center, histHOA_ver2_DipX_Asym))


histFN_DipY_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=Amplitude_1_Asym/normalization_01_Asym)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleY-Asym.hist',zip(bin_center, histFN_DipY_Asym))

histHOA_DipY_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleY-Asym.hist',zip(bin_center, histHOA_DipY_Asym))

histHOA_ver2_DipY_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=Asym_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(Asym_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleY-Asym-ver2.hist',zip(bin_center, histHOA_ver2_DipY_Asym))


histFN_DipZ_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=Amplitude_1_Asym/normalization_01_Asym)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleZ-Asym.hist',zip(bin_center, histFN_DipZ_Asym))

histHOA_DipZ_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleZ-Asym.hist',zip(bin_center, histHOA_DipZ_Asym))

histHOA_ver2_DipZ_Asym,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=Asym_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(Asym_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleZ-Asym-ver2.hist',zip(bin_center, histHOA_ver2_DipZ_Asym))



#---ZDisp---

histFN_DipX_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=Amplitude_1_ZDisp/normalization_01_ZDisp)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleX-ZDisp.hist',zip(bin_center, histFN_DipX_ZDisp))

histHOA_DipX_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleX-ZDisp.hist',zip(bin_center, histHOA_DipX_ZDisp))

histHOA_ver2_DipX_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False, weights=ZDisp_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(ZDisp_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleY-ZDisp-ver2.hist',zip(bin_center, histHOA_ver2_DipX_ZDisp))

histFN_DipY_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=Amplitude_1_ZDisp/normalization_01_ZDisp)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleY-ZDisp.hist',zip(bin_center, histFN_DipY_ZDisp))

histHOA_DipY_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleY-ZDisp.hist',zip(bin_center, histHOA_DipY_ZDisp))

histHOA_ver2_DipY_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False, weights=ZDisp_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(ZDisp_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleY-ZDisp-ver2.hist',zip(bin_center, histHOA_ver2_DipY_ZDisp))

histFN_DipZ_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=Amplitude_1_ZDisp/normalization_01_ZDisp)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('FN-'+fileNameString+'-DipoleZ-ZDisp.hist',zip(bin_center, histFN_DipZ_ZDisp))

histHOA_DipZ_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleZ-ZDisp.hist',zip(bin_center, histHOA_DipZ_ZDisp))

histHOA_ver2_DipZ_ZDisp,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False, weights=ZDisp_coord[nonZeros_0]*DW0_0[nonZeros_0]/(np.sqrt(np.sum(ZDisp_coord[nonZeros_0]**2*DW0_0[nonZeros_0]))*normQ00))
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('HOA-'+fileNameString+'DipoleZ-ZDisp-ver2.hist',zip(bin_center, histHOA_ver2_DipZ_ZDisp))




histDipoleY,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,1],bins=NBins, normed=False,weights=DW0_0[nonZeros_0]/normQ00**2)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('DipoleY.hist',zip(bin_center,histDipoleY))

histDipoleX,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,0],bins=NBins, normed=False,weights=DW0_0[nonZeros_0]/normQ00**2)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('DipoleX.hist',zip(bin_center,histDipoleX))

histDipoleZ,bin_edges=np.histogram(DipoleMoment_0[nonZeros_0][:,2],bins=NBins, normed=False,weights=DW0_0[nonZeros_0]/normQ00**2)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
np.savetxt('DipoleZ.hist',zip(bin_center,histDipoleZ))

np.savetxt('muxvsq.data',zip(DipoleMoment_0[nonZeros_0][:,0],Asym_coord[nonZeros_0],ZDisp_coord[nonZeros_0],Q_Asym,Q_ZDisp))
np.savetxt('muxvsAmp.data',zip(DipoleMoment_0[nonZeros_0][:,0],Amplitude_1_Asym/normalization_01_Asym,Amplitude_1_ZDisp/normalization_01_ZDisp,Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00),Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00)))

np.savetxt('muyvsq.data',zip(DipoleMoment_0[nonZeros_0][:,1],Asym_coord[nonZeros_0],ZDisp_coord[nonZeros_0],Q_Asym,Q_ZDisp))
np.savetxt('muyvsAmp.data',zip(DipoleMoment_0[nonZeros_0][:,1],Amplitude_1_Asym/normalization_01_Asym,Amplitude_1_ZDisp/normalization_01_ZDisp,Q_Asym*DW0_0[nonZeros_0]/(normQAsym*normQ00),Q_ZDisp*DW0_0[nonZeros_0]/(normQZDisp*normQ00)))
