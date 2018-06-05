#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys
import glob

def calcMass(self,particle,name):
    dtau=self.Wfn.dtau

    massH=1.00782503223
    massD=2.0141017778
    if particle=='H':
        massHamu=massH*massConversionFactor
    elif particle=='D':
        massHamu=massD*massConversionFactor
    if name=='SharedProton':
        U=x[:,1,:]-x[:,0,:]
        V=x[:,3,:]-x[:,0,:]
        magU=np.sqrt(U[:,0]**2+U[:,1]**2+U[:,2]**2)
        magV=np.sqrt(V[:,0]**2+V[:,1]**2+V[:,2]**2)
        costheta= np.diag(np.dot(U,V.T))/(magU*magV)
        mass=1.0/(2.0*((1.0/(massOamu))+((1.000000-costheta)/(massHamu))))
        
        dist=self.Wfn.molecule.function(x)
        
    



    return mass



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

plt.figure(2)
dtau=10.0
dist=R_dr
Wfn.setNodalSurface('SharedProton')
Wfn.setIsotope('notDeuterated')
mass=Wfn.molecule.calcReducedmass(coords)
P_rec=np.exp(-2.0*dist*dist*mass/dtau)
plt.scatter(R_dr,P_rec,color='blue')



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

plt.figure(2)
dtau=10.0
dist=R_dr
Wfn.setNodalSurface('SharedProton')
Wfn.setIsotope('fullyDeuterated')
mass=Wfn.molecule.calcReducedmass(coords)
P_rec=np.exp(-2.0*dist*dist*mass/dtau)
plt.scatter(R_dr,P_rec,color='green')





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


plt.figure(2)
dtau=10.0
dist=R_dr
Wfn.setNodalSurface('Z-displacement')
Wfn.setIsotope('notDeuterated')
mass=Wfn.molecule.calcReducedmass(coords)
P_rec=np.exp(-2.0*dist*dist*mass/dtau)
plt.scatter(R_dr,P_rec,color='teal')




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

plt.figure(2)
dtau=10.0
dist=R_dr
Wfn.setNodalSurface('Z-displacement')
Wfn.setIsotope('fullyDeuterated')
mass=Wfn.molecule.calcReducedmass(coords)
P_rec=np.exp(-2.0*dist*dist*mass/dtau)
plt.scatter(R_dr,P_rec,color='blue')



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


plt.show()
