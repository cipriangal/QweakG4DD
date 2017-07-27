#!/apps/python/PRO/bin/python
from subprocess import call
import sys,os,time

def main():
    
    directory="/lustre/expphy/volatile/hallc/qweak/ciprian/farmoutput/g41001p01/sample/moustaches/withShower/md3"
    nEv=100000
    nrStop=600
    nrStart=0
    pol="V"
    nDist=203

    idRoot = pol+'_sampled'+str(nDist)+'_%03dk'% (nEv/1000)

    for nr in range(nrStart,nrStop): # repeat for nr jobs
        idN = idRoot+'_%05d'% (nr)
        print idN
        
        if not os.path.exists(directory+"/"+idN+"/log"):
            os.makedirs(directory+"/"+idN+"/log")

        if pol=="V":
            call("root -l -q -b ../rootScripts/samplePrimaryDist.C\\("+str(nEv)+",1,"+str(nr)+","+str(nDist)+"\\)",shell=True)
        else:
            call("root -l -q -b ../rootScripts/samplePrimaryDist.C\\("+str(nEv)+",-1,"+str(nr)+","+str(nDist)+"\\)",shell=True)
        
        call(["mv","positionMomentum_%d.in"%(nr),directory+"/"+idN+"/positionMomentum.in"])
        call(["mv","polarization_%d.in"%(nr),directory+"/"+idN+"/polarization.in"])

        
    print "I am all done"

                    
if __name__ == '__main__':
    main()
                            
