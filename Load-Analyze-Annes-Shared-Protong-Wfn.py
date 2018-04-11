#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import glob

if len(sys.argv)>1:
    print 'Usage: just use it, no arguments'
    end



starttime=time.time()
au2wn=219474.63
nBins=51
ang2bohr=1.88973
min,max=0,2.5


Wfn=dmc.wavefunction('HOHOH', 20000)
Destination='AnnesH3O2/'

coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_es_g.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)
print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val


R_dr=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_zdisp=Wfn.molecule.calcZdisplacement(coords)
#R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
#R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)



HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)

plt.subplot(211)
plt.plot(bin_center,HistDR,color='green',linewidth=4,linestyle='-')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='green',linewidth=4,linestyle='-')


coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_es_g-alld.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)
print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val


R_dr=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_zdisp=Wfn.molecule.calcZdisplacement(coords)
#R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
#R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)



HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='green',linewidth=4,linestyle='--')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='green',linewidth=4,linestyle='--')





coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_ez.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)
print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val


R_dr=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_zdisp=Wfn.molecule.calcZdisplacement(coords)
#R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
#R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='blue',linewidth=4,linestyle='-')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='blue',linewidth=4,linestyle='-')


coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_ez-alld.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)
print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val


R_dr=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_zdisp=Wfn.molecule.calcZdisplacement(coords)
#R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
#R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='blue',linewidth=4,linestyle='--')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='blue',linewidth=4,linestyle='--')





fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedZ-Displacement-2000-5-25-5Eq-*')
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
print ''
R_dr=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
R_zdisp=Wfn.molecule.calcZdisplacement(gathercoords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='red',linewidth=4,linestyle='-')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='red',linewidth=4,linestyle='-')




fileList=glob.glob('Wfn-DODOD-Tau/DODOD-Start-Ground-Propagate-ZDisp-DW-ZDisp2000*')
print '       hopefully the z displacement, deuterated??', fileList
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
print ''
R_dr=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
R_zdisp=Wfn.molecule.calcZdisplacement(gathercoords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='red',linewidth=4,linestyle='--')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='red',linewidth=4,linestyle='--')






#fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-2000-5-50*')


fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-Excited-DeltaR-2000-5-25-5Eq-*')
#fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-2000-5-25-5Eq-*')
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
print ''
R_dr=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
R_zdisp=Wfn.molecule.calcZdisplacement(gathercoords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='magenta',linewidth=4,linestyle='-')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='magenta',linewidth=4,linestyle='-')


#fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-2000-5-50*')

fileList=glob.glob("Wfn-DODOD-Tau/DODOD-Start-Ground-Propagate-SharedProton-DW-SharedProton2000-*")
#fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-2000-5-25-5Eq-*')
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
print ''
R_dr=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
R_zdisp=Wfn.molecule.calcZdisplacement(gathercoords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.plot(bin_center,HistDR,color='magenta',linewidth=4,linestyle='--')
plt.subplot(212)
plt.plot(bin_center,HistZD,color='magenta',linewidth=4,linestyle='--')




fileList=glob.glob("Wfn-DODOD-Tau/DODOD-Start-Ground-Propagate-SharedProton-DW-SharedProton20000-*")
#fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-ExcitedRn-2000-5-25-5Eq-*')
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
print ''
R_dr=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
R_zdisp=Wfn.molecule.calcZdisplacement(gathercoords)

HistDR,bin_edges_rock=np.histogram(np.absolute(R_dr),bins=nBins,range=(min,max),density=True)
HistZD,bin_edges_rock=np.histogram(np.absolute(R_zdisp),bins=nBins,range=(min,max),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(211)
plt.xlim(0.000, 0.15)
plt.ylim(0.000,0.4)
plt.plot(bin_center,HistDR,color='black',linewidth=4,linestyle='--')
plt.subplot(212)
plt.xlim(0.000, 0.15)
plt.ylim(0.0000,.6)
plt.plot(bin_center,HistZD,color='black',linewidth=4,linestyle='--')



plt.xlim(0.000, 0.25)

plt.show()
