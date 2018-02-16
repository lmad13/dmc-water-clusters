import DMCClusters as dmc
Wfn=dmc.wavefunction('HOHOH', 5)
print 'Before: Wfn.x',Wfn.x
Wfn.propagate(Wfn.x,5)
print 'After: Wfn.x',Wfn.x
print 'done!'
