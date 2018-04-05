#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
import sys


if len(sys.argv)>1:
    print 'Usage: just use it, no arguments'
    end



starttime=time.time()
au2wn=219474.63
nBins=51
ang2bohr=1.88973

asMin,asMax=0,0.75
rspMin,rspMax=0,2.0
rockMin,rockMax=1.0,5.0

Wfn=dmc.wavefunction('HOHOH', 20000)
Destination='AnnesH3O2/'

coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_e.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)
print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val

R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
R_sp=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)
max_is_here=np.argmax(R_as)

print 'largest', R_as[max_is_here],'is here',max_is_here,'with coords \n', coords[max_is_here]

print swapped

HistRock1,bin_edges_rock=np.histogram(R_rock1,bins=nBins,range=(rockMin,rockMax),density=True)
HistRock3,bin_edges_rock=np.histogram(R_rock3,bins=nBins,range=(rockMin,rockMax),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(313)
plt.plot(bin_center,HistRock1,color='red')
plt.plot(bin_center,HistRock3,color='magenta')


HistTemp,bin_edges=np.histogram(np.absolute(R_as),bins=nBins,range=(asMin,asMax),density=True)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.subplot(312)
plt.plot(bin_center,HistTemp,linewidth=4,color='red')#,color=linecolor[iwfn])  

plt.subplot(311)
HistTempRsp,bin_edges_rsp=np.histogram(np.absolute(R_sp),bins=nBins,range=(rspMin,rspMax),density=True)
bin_center_rsp=(bin_edges_rsp[:-1]+bin_edges_rsp[1:])/2.0
plt.plot(bin_center_rsp,HistTempRsp,linewidth=4,color='red')#,color=linecolor[iwfn])  

Histas_rsp,xedges_as,yedges_rsp=np.histogram2d(R_as,R_sp,bins=nBins,range=[[asMin,asMax],[rspMin,rspMax]],normed=True)
plt.figure(2)
plt.matshow(np.log(Histas_rsp))

Histas_rock,xedges_as,yedges_rock=np.histogram2d(R_as,R_rock1,bins=nBins,range=[[asMin,asMax],[rockMin,rockMax]],normed=True)
plt.figure(5)
plt.matshow(np.log(Histas_rock))



###
###
###

coords,dw=Wfn.loadCoords('AnnesH3O2/walkers_e2.xyz')
print 'shape of coords', coords.shape
coords=coords/ang2bohr
print 'average dw',np.average(dw)
InternalNames,Values=Wfn.molecule.calcAverageInternalCoordinates(coords,dw)

print 'Internal Coordinates! Min, Max, avwerage, expecatation value, std'
for iname,val in zip(InternalNames,Values):
    print iname, val

R_as,swapped=Wfn.molecule.calcStretchAnti(coords)
R_sp=Wfn.molecule.calcSharedProtonDisplacement(coords)
R_rock1,R_rock3=Wfn.molecule.calcRocks(coords)
max_is_here=np.argmax(R_as)
print 'largest', R_as[max_is_here],'is here',max_is_here,'with coords \n', coords[max_is_here]

print swapped


HistRock1,bin_edges_rock=np.histogram(R_rock1,bins=nBins,range=(rockMin,rockMax),density=True)
HistRock3,bin_edges_rock=np.histogram(R_rock3,bins=nBins,range=(rockMin,rockMax),density=True)
bin_center=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0
plt.figure(1)
plt.subplot(313)
plt.plot(bin_center,HistRock1,color='blue')
plt.plot(bin_center,HistRock3,color='#0504aa')


HistTemp,bin_edges=np.histogram(np.absolute(R_as),bins=nBins,range=(asMin,asMax),density=True)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.subplot(312)
plt.plot(bin_center,HistTemp,linewidth=4,color='blue')#,color=linecolor[iwfn])                                                                                          
plt.subplot(311)
HistTempRsp,bin_edges_rsp=np.histogram(np.absolute(R_sp),bins=nBins,range=(rspMin,rspMax),density=True)
bin_center_rsp=(bin_edges_rsp[:-1]+bin_edges_rsp[1:])/2.0
plt.plot(bin_center_rsp,HistTempRsp,linewidth=4,color='blue')#,color=linecolor[iwfn])      

Histas_rsp,xedges_as,yedges_rsp=np.histogram2d(R_as,R_sp,bins=nBins,range=[[asMin,asMax],[rspMin,rspMax]],normed=True)
plt.figure(3)
plt.matshow(np.log(Histas_rsp))

Histas_rock,xedges_as,yedges_rock=np.histogram2d(R_as,R_rock1,bins=nBins,range=[[asMin,asMax],[rockMin,rockMax]],normed=True)
plt.figure(6)
plt.matshow(np.log(Histas_rock))

###
###
###

Hist=np.zeros((nBins))
HistRsp=np.zeros((nBins))
HistRock1=np.zeros((nBins))
HistRock3=np.zeros((nBins))
gathercoords=[]
import glob
fileList=glob.glob('Wfn-HOHOH-Tau/HOHOH-Start-Ground-Propagate-Asym-DW-Asym-20000-5-25-*')
for i,fname in enumerate(fileList):
    print 'loading', fname,
    coords,dw=Wfn.loadCoords(fname)
    if i==0:
        gathercoords=coords
    else:
        gathercoords=np.concatenate((gathercoords,coords),axis=0)
R_as,swapped=Wfn.molecule.calcStretchAnti(gathercoords)
R_sp=Wfn.molecule.calcSharedProtonDisplacement(gathercoords)
                                      
R_rock1,R_rock3=Wfn.molecule.calcRocks(gathercoords)
Hist,bin_edges=np.histogram(np.absolute(R_as),bins=nBins,range=(asMin,asMax),density=True)
bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0


HistRsp,bin_edges_rsp=np.histogram(np.absolute(R_sp),bins=nBins,range=(rspMin,rspMax),density=True)
bin_center_rsp=(bin_edges_rsp[:-1]+bin_edges_rsp[1:])/2.0

HistRock1,bin_edges_rock=np.histogram(R_rock1,bins=nBins,range=(rockMin,rockMax),density=True)
HistRock3,bin_edges_rock=np.histogram(R_rock3,bins=nBins,range=(rockMin,rockMax),density=True)
bin_center_rock=(bin_edges_rock[:-1]+bin_edges_rock[1:])/2.0

plt.figure(1)
plt.subplot(313)
plt.plot(bin_center_rock,HistRock1,color='green')
plt.plot(bin_center_rock,HistRock3,color='#25a36f')


plt.subplot(311)
plt.plot(bin_center_rsp,HistRsp,linewidth=4,color='green')#,color=linecolor[iwfn])      
plt.subplot(312)
plt.plot(bin_center,Hist,linewidth=4,color='green')#,color=linecolor[iwfn])                     


plt.figure(4)
Histas_rsp,xedges_as,yedges_rsp=np.histogram2d(R_as,R_sp,bins=nBins,range=[[asMin,asMax],[rspMin,rspMax]],normed=True)
plt.matshow(np.log(Histas_rsp))

Histas_rock,xedges_as,yedges_rock=np.histogram2d(R_as,R_rock1,bins=nBins,range=[[asMin,asMax],[rockMin,rockMax]],normed=True)
plt.figure(6)
plt.matshow(np.log(Histas_rock))

endtime=time.time()
print 'that took', (endtime-starttime)/60.0 , 'min'

plt.show()
plt.savefig('AnnesH3O2/histogram.png')



print 'done!'
