#!/apps/python/PRO/bin/python
from subprocess import call
import sys,os,time

def main():

    #beamE=[5,10,30,50,100]
    beamEnergy=[10,50]

    #xPos=[-0.1,-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05,0.1]
    xPos=[0,60]

    xAng=[-40,-20,-10,-5,0]
    #xAng=[-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,-1,0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]

    yPos=[326,332,335,338,344]

    #yAng=[-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,-1,0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]
    yAng=[-80,-40,-20,-10,0,5,15,30,60]

    pol="L"

    email="ciprian@jlab.org"
    source="/lustre/expphy/work/hallc/qweak/ciprian/simCodeG410/QweakG4DD"
    directory="/lustre/expphy/volatile/hallc/qweak/ciprian/farmoutput/acrossYposCheck"

    nEv=20000
    nrStart = 0
    nrStop  = 1
    submit=1

    for xA in xAng: # angle along x axis (along MD3, perp MD1)
        for xP in xPos: # x position of the beam
            for beamE in beamEnergy: # E of the beam
                for nr in range(nrStart,nrStop): # repeat for nr jobs
                    for yP in yPos:
                        for yA in yAng:
                            zP=575.0 #in front of MD no Pb (this should go to at most 88 deg until the end of the bar)
                            idN= pol+'_%04d_%06.2f_%06.2f_%06.2f_%06.2f_%06.2f_%03d'% (beamE,xP,yP,zP,xA,yA,nr)
                            createMacFile(directory,idN,xP,yP,zP,xA,yA,beamE,pol,nEv,nr)
                            createXMLfile(idN,directory,email,source)
                            call(["cp",source+"/build/QweakSimG4",directory+"/"+idN+"/QweakSimG4"])
                            call(["cp",source+"/build/QweakGetPEs",directory+"/"+idN+"/QweakGetPEs"])
                            call(["cp",source+"/myQweakCerenkovOnly.mac",directory+"/"+idN+"/myQweakCerenkovOnly.mac"])
                            #sys.exit()#for testing purposes

                            if submit==1:
                                print "submitting X position", xP,"with angle",xA,"for the ",nr," time"," with E ",beamE
                                call(["jsub","-xml",directory+"/"+idN+"/job.xml"])
                            else:
                                print "do not submit ",submit


def createMacFile(directory,idname,xPos,yPos,zPos,xAng,yAng,beamE,pol,nEv,nr):
    if not os.path.exists(directory+"/"+idname):
        os.makedirs(directory+"/"+idname)

    f=open(directory+"/"+idname+"/myRun.mac",'w')
    f.write("/control/execute myQweakCerenkovOnly.mac\n")
    f.write("/PrimaryEvent/SetBeamPositionX "+str(xPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionY "+str(yPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionZ "+str(zPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamDirectionX "+str(xAng)+" deg\n")
    f.write("/PrimaryEvent/SetBeamDirectionY "+str(yAng)+" deg\n")
    f.write("/PrimaryEvent/SetFixedPosMom true\n")
    f.write("/PrimaryEvent/SetPolarization "+pol+"\n")
    f.write("/EventGen/SetBeamEnergy    "+str(beamE)+" MeV\n")
    f.write("/TrackingAction/TrackingFlag 3\n")
    f.write("/EventGen/SelectOctant 3\n")
    seedA=int(time.time()/2000.) +   100*nr + nr
    seedB=int(time.time()/300.)  + 10000*nr + nr
    f.write("/random/setSeeds "+str(seedA)+" "+str(seedB)+"\n")
    f.write("/run/beamOn "+str(nEv)+"\n")
    f.close()
    return 0

def createXMLfile(idname,directory,email,source):
    if not os.path.exists(directory+"/"+idname+"/log"):
        os.makedirs(directory+"/"+idname+"/log")

    f=open(directory+"/"+idname+"/job.xml","w")
    f.write("<Request>\n")
    f.write("  <Email email=\""+email+"\" request=\"false\" job=\"true\"/>\n")
    f.write("  <Project name=\"qweak\"/>\n")
    f.write("  <Track name=\"simulation\"/>\n")
    f.write("  <Name name=\""+idname+"\"/>\n")
    f.write("  <OS name=\"centos7\"/>\n")
    f.write("  <Command><![CDATA[\n")
    f.write("cd "+directory+"/"+idname+"\n")
    f.write("./QweakSimG4 myRun.mac\n")
    f.write("./QweakGetPEs 3 QwSim_0.root\n")
    f.write("  ]]></Command>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")
    #f.write("  <TimeLimit unit=\"minutes\" time=\"4320\"/>\n")
    f.write("  <Job>\n")
    f.write("    <Stdout dest=\""+directory+"/"+idname+"/log/log.out\"/>\n")
    f.write("    <Stderr dest=\""+directory+"/"+idname+"/log/log.err\"/>\n")
    f.write("  </Job>\n")
    f.write("</Request>\n")
    f.close()

    return 0

if __name__ == '__main__':
    main()

