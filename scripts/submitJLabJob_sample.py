#!/apps/python/PRO/bin/python
from subprocess import call
import sys,os,time

def main():
    
    #center, x,y,z=0,335,560
    _xP=0.
    _yP=335.0
    _zP=560.
    _Px=0.#deg
    _Py=0.
    _beamE=1160#MeV
    _tracking=0#0=primary only | 1=prim + opt photon | 2=no optical ph and 10x faster than 3=full
    _email="ciprian@jlab.org"
    _source="/lustre/expphy/work/hallc/qweak/ciprian/simCodeG410/QweakG4DD"
    _directory="/lustre/expphy/volatile/hallc/qweak/ciprian/farmoutput/g41001p01/sampled/mc3MtCalc_0001mm_v2"
    _nEv=5000
    _nrStop=1000
    _nrStart=0
    _pol="V"
    submit=1

    idRoot= _pol+'_sampled_%03dk'% (_nEv/1000) 
    for nr in range(_nrStart,_nrStop): # repeat for nr jobs
        _idN= idRoot+'_%04d'% (nr) 
        print _idN
        createMacFile(_directory,_idN,_xP,_yP,_zP,_Px,_Py,_tracking,_beamE,_pol,_nEv,nr)
        ##create input files
        seedA=int(time.time()/1346.)+10000*nr+nr
        if _pol=="V":
            call("root -l -q -b ../rootScripts/samplePrimaryDist.C\\("+str(seedA)+","+str(_nEv)+",1\\)",shell=True)
        else:
            call("root -l -q -b ../rootScripts/samplePrimaryDist.C\\("+str(seedA)+","+str(_nEv)+",-1\\)",shell=True)
        call(["mv","positionMomentum.in",_directory+"/"+_idN+"/positionMomentum.in"])
        call(["mv","polarization.in",_directory+"/"+_idN+"/polarization.in"])

        call(["cp",_source+"/build/QweakSimG4",_directory+"/"+_idN+"/QweakSimG4"])
        call(["cp",_source+"/myQweakCerenkovOnly.mac",_directory+"/"+_idN+"/myQweakCerenkovOnly.mac"])


    createXMLfile(_source,_directory,idRoot,_nrStart,_nrStop,_email)


    if submit==1:
        print "submitting position sampled with id",_idN," between ",_nrStart,_nrStop
        call(["jsub","-xml",_source+"/scripts/jobs/"+idRoot+".xml"])
    else:
        print "NOT submitting position sampled with id",_idN," between ",_nrStart,_nrStop
        
    print "I am all done"


def createMacFile(directory,idname,
                  xPos,yPos,zPos,
                  Px,Py,tracking,
                  beamE,pol,nEv,nr):
    if not os.path.exists(directory+"/"+idname+"/log"):
        os.makedirs(directory+"/"+idname+"/log")
   
    f=open(directory+"/"+idname+"/myRun.mac",'w')
    f.write("/control/execute myQweakCerenkovOnly.mac\n")
    f.write("/PrimaryEvent/SetBeamPositionX "+str(xPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionY "+str(yPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionZ "+str(zPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamDirectionX "+str(Px)+" deg\n")
    f.write("/PrimaryEvent/SetBeamDirectionY "+str(Py)+" deg\n")
    f.write("/PrimaryEvent/SetFixedPosMom false\n")
    f.write("/PrimaryEvent/SetPolarization "+pol+"\n")
    f.write("/EventGen/SetBeamEnergy    "+str(beamE)+" MeV\n")
    f.write("/TrackingAction/TrackingFlag "+str(tracking)+"\n")
    f.write("/EventGen/SelectOctant 3\n")
    seedA=int(time.time()/2000.)+   100*nr+nr
    seedB=int(time.time()/300. ) +10000*nr+nr
    f.write("/random/setSeeds "+str(seedA)+" "+str(seedB)+"\n")
    f.write("/run/beamOn "+str(nEv)+"\n")
    f.close()
    return 0

def createXMLfile(source,writeDir,idRoot,nStart,nStop,email):
    
    if not os.path.exists(source+"/scripts/jobs"):
        os.makedirs(source+"/scripts/jobs")

    f=open(source+"/scripts/jobs/"+idRoot+".xml","w")
    f.write("<Request>\n")
    f.write("  <Email email=\""+email+"\" request=\"false\" job=\"true\"/>\n")
    f.write("  <Project name=\"qweak\"/>\n")
#    f.write("  <Track name=\"debug\"/>\n")
    f.write("  <Track name=\"simulation\"/>\n")
    f.write("  <Name name=\""+idRoot+"\"/>\n")
    f.write("  <OS name=\"centos65\"/>\n")
    f.write("  <Command><![CDATA[\n")
    f.write("QweakSimG4 myRun.mac\n")
    f.write("  ]]></Command>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")

    for nr in range(nStart,nStop): # repeat for nr jobs
        f.write("  <Job>\n")
        idName= writeDir+"/"+idRoot+'_%04d'%(nr)
        f.write("    <Input src=\""+idName+"/QweakSimG4\" dest=\"QweakSimG4\"/>\n")
        f.write("    <Input src=\""+idName+"/myQweakCerenkovOnly.mac\" dest=\"myQweakCerenkovOnly.mac\"/>\n")
        f.write("    <Input src=\""+idName+"/myRun.mac\" dest=\"myRun.mac\"/>\n")
        f.write("    <Input src=\""+idName+"/positionMomentum.in\" dest=\"positionMomentum.in\"/>\n")
        f.write("    <Input src=\""+idName+"/polarization.in\" dest=\"polarization.in\"/>\n")
        f.write("    <Output src=\"QwSim_0.root\" dest=\""+idName+"/QwSim_0.root\"/>\n")
        f.write("    <Output src=\"o_tuple.root\" dest=\""+idName+"/o_tuple.root\"/>\n")
        f.write("    <Stdout dest=\""+idName+"/log/log.out\"/>\n")
        f.write("    <Stderr dest=\""+idName+"/log/log.err\"/>\n")
        f.write("  </Job>\n\n")

    f.write("</Request>\n")
    f.close()
    return 0

                    
if __name__ == '__main__':
    main()
                            
