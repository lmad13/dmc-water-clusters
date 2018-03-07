import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys

print 'testing HOHOH-'

N_size=10
propagationSteps=5
au2wn=219474.63

Wfn=dmc.wavefunction('HOHOH', N_size)

initialx=Wfn.x*1.01

print 'testing initial V ............',np.all(np.isclose(Wfn.molecule.V(initialx)*au2wn,540.54408768))

v_ref_equilibration,pop_list_equilibration,equilibratedCoordinates,descendants=Wfn.propagate(initialx,propagationSteps,printCensus=True,initialPop=N_size,testing=True)
print 'propagating ..................'
idealV_ref=[0.002462900097768067, 0.007253169644062676, 0.010156689336005748, 0.012494787956209688, 0.017783383299700296, 0.02318534702796856]
print 'testing vref .................', np.all(np.isclose(v_ref_equilibration,idealV_ref))

idealRn=[-0.07655971,  0.27293458, -0.24240796,  0.11960364, -0.17116614,  0.09562018, 0.30075851, -0.02536788,  0.05155765]
print 'testing Rn ...................',np.all(np.isclose(Wfn.molecule.calcSharedProtonDisplacement(equilibratedCoordinates),idealRn,rtol=1e-6))


print 'testing Eckart Rotation with Rn ............',

eckRotx=Wfn.molecule.eckartRotate(equilibratedCoordinates)

print np.all(np.isclose(Wfn.molecule.calcSharedProtonDisplacement(eckRotx),idealRn, rtol=1e-6))

idealDipole=[[ 0.19743043 , 0.92138308, -0.11961317],
             [-0.28923288 , 0.87957105,  0.10387864],
             [ 0.17616207 , 0.96483681, -0.05132077],
             [-0.0826187  , 0.94220673, -0.01086911],
             [ 0.18083071 , 0.94145279,  0.03618015],
             [-0.03283231 , 0.91121847, -0.15905377],
             [-0.2900713  , 0.95461155, -0.11632933],
             [-0.0024041  , 0.97352727,  0.11038614],
             [-0.00320079 , 0.90135676, -0.06295138]]

print 'testing calcDipole .........................',np.all(np.isclose(Wfn.molecule.calcDipole(equilibratedCoordinates),idealDipole,rtol=1e-6))

Wfn.exportCoords(equilibratedCoordinates,'testingEquilibratedCoordinates.xyz',descendants)
Wfn.exportCoords(eckRotx,'testingEckartRotatedCoordinates.xyz',descendants)
