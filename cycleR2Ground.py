import sys
import os
import time
starttime=time.time()
count=0
timeRate=1.11e-06
nSteps=500
for nDescSteps in [25,50,100,150]:#range(25,75,25):
    for nRepsDesc in range(5,20,10):
        for nReps in range (5,20,10):
            for nWalkers in [20000,10000,5000,20000]:#range(2000,3000,2000):

                print '\n\n################    WORKING ON #',count,
                count=count+1
                print 'nWalkers:   ',nWalkers
                print 'nReps:      ',nReps
                print 'nRepsDesc:  ',nRepsDesc
                print 'nDescSteps:  ',nDescSteps
                print 'nSteps:     ',nSteps
                print '#####   Estimated Time:',str(timeRate*nWalkers*nReps*nDescSteps*nRepsDesc*nSteps),'s or ',str(timeRate*nWalkers*nReps*nDescSteps*nRepsDesc*nSteps/60.0),'min #####\n'
                executeStartTime=time.time()
                sys.stdout.flush()
                os.system('./HOHOH-groundState.py '+str(nWalkers)+' '+str(nReps)+
                          ' '+str(nDescSteps)+' '+str(nRepsDesc))
                executionTime=time.time()-executeStartTime
                timeRate=executionTime/(nWalkers*nReps*nDescSteps*nRepsDesc*nSteps)
                print 'that took', executionTime, 's with a timerate of ', timeRate
                sys.stdout.flush()

print 'done cycling!', time.time()-starttime

