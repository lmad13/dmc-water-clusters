import numpy as np
import matplotlib.pyplot as plt
import DMCClusters as dmc
import time
starttime=time.time()
au2wn=219474.63
averaged_vref=[]
list_of_pop_list=[]
N_size=5000
N_steps=300
nReps=10
Wfn=dmc.wavefunction('HOHOH', N_size)
initialx=Wfn.x
print 'initial V', Wfn.molecule.V([initialx[0]])*au2wn

for iwfn in range(nReps):
    v_ref_list,pop_list,finalCoords,descendants=Wfn.propagate(initialx,N_steps,printCensus=False,initialPop=N_size)
    endtime=time.time()
    averaged_vref.append(np.average(np.array(v_ref_list)[N_steps-100:N_steps])*au2wn)
    list_of_pop_list.append(pop_list)
    plt.plot(np.array(v_ref_list)*au2wn)
    initialx=finalCoords
plt.show()
print 'averaged v_ref:',averaged_vref
print 'the average of average V_ref for the last 100 steps is',np.average(np.array(averaged_vref)), ' cm-1',
print 'standard deviation', np.std(np.array(averaged_vref)), ' cm-1'
print 'uncertainity is', (np.max(averaged_vref)-np.min(averaged_vref))/(2.0*np.sqrt(nReps))
for pop_list in list_of_pop_list:
    plt.plot(pop_list)
plt.show()



print 'that took', endtime-starttime, 'seconds and ', (endtime-starttime)/60.0 , 'minutes'


print 'done!'
