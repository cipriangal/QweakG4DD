#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "QweakSimUserMainEvent.hh"
#include "QweakSimSystemOfUnits.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"

using namespace std;

int evNr;
int primary;//0 secondary, 1 primary
float x,y,z,E;
float angX,angY,angXi,angYi;
float xi,yi,zi;
float polT,polTi;
double asymPpM,asymPmM;
std::vector<double> xI,yI,zI,aXi,aYi,polI;

void findInt(std::vector<int> &inter,std::vector<int> &val, int trackID,int parent, int &hasPar, int &nInt);
void processOne(TTree *QweakSimG4_Tree,TTree *tout, int &nrEvts);
void readInitial(string fnm);


int main(int argc, char** argv){

  if( argc < 2 ) {
    cout<<" usage: build/QweakAna [path to infile with list of output QweakSimG4 trees]"<<endl;
    return 1;
  }

  string files(argv[1]);

  TFile *fout=new TFile("o_hits.root","RECREATE");
  TH1I *hEntries=new TH1I("hEntries","number of processed events",2,0,2);
  TTree *tout=new TTree("t","Stripped QweakSimG4 tree for hits");
  tout->Branch("evNr",&evNr,"evNr/I");
  tout->Branch("primary",&primary,"primary/I");
  tout->Branch("x",&x,"x/F");
  tout->Branch("y",&y,"y/F");
  tout->Branch("z",&z,"z/F");
  tout->Branch("E",&E,"E/F");
  tout->Branch("angX",&angX,"angX/F");
  tout->Branch("angY",&angY,"angY/F");
  tout->Branch("polT",&polT,"polT/F");
  tout->Branch("xi",&xi,"xi/F");
  tout->Branch("yi",&yi,"yi/F");
  tout->Branch("zi",&zi,"zi/F");
  tout->Branch("angXi",&angXi,"angXi/F");
  tout->Branch("angYi",&angYi,"angYi/F");
  tout->Branch("polTi",&polTi,"polTi/F");
  tout->Branch("asymPpM",&asymPpM,"asymPpM/D");
  tout->Branch("asymPmM",&asymPmM,"asymPmM/D");
  
  int totEv=0;
  int simEvts(0);
  if ( files.find(".root") < files.size() ){
    readInitial(files);
    TFile *fin=new TFile(files.c_str(),"READ");
    if(!fin->IsOpen()){
      cout<<"Problem: can't find file: "<<files<<endl;
      fout->Close();
      return -3;
    }

    TTree *QweakSimG4_Tree=(TTree*)fin->Get("QweakSimG4_Tree");
    cout<<"processing only one file: "<<files.c_str()<<endl;    
    
    int evts=QweakSimG4_Tree->GetEntries();
    totEv+=evts;
    cout<<" total nr ev: "<<evts<<" "<<totEv<<endl;
    
    processOne(QweakSimG4_Tree,tout,simEvts);
    simEvts=std::ceil(simEvts/1000.)*1000;
    
    cout<<"Done looping!"<<endl;
    fin->Close();
    delete fin;
  }else{
    ifstream ifile(files.c_str());
    string data;
    while(ifile>>data){
      readInitial(data);
      TFile *fin=new TFile(data.c_str(),"READ");
      if(!fin->IsOpen()){
	cout<<"Problem: can't find file: "<<data<<endl;
	continue;
      }
      TTree *QweakSimG4_Tree=(TTree*)fin->Get("QweakSimG4_Tree");
      cout<<"processing : "<<data.c_str()<<endl;          
      int evts=QweakSimG4_Tree->GetEntries();
      totEv+=evts;
      cout<<" total nr ev: "<<evts<<" "<<totEv<<endl;      
      processOne(QweakSimG4_Tree,tout,simEvts);      
      simEvts=std::ceil(simEvts/1000.)*1000;
      fin->Close();
      delete fin;
    }
    ifile.close();
  }
  
  cout<<"Processed "<<totEv<<" events"<<endl;
  cout<<"\tfrom "<<simEvts<<" simulated events"<<endl;
  
  hEntries->SetBinContent(1,totEv);
  hEntries->SetBinContent(2,simEvts);
  fout->cd();
  hEntries->Write();
  tout->Write();
  fout->Close();
  return 0;
}

void processOne(TTree *QweakSimG4_Tree, TTree *tout, int &nrEvts){

  std::vector<int> interaction;
  std::vector<int> trackID;

  //set addresses of leafs
  QweakSimUserMainEvent* event = 0;
  QweakSimG4_Tree->SetBranchAddress("QweakSimUserMainEvent",&event);    

  for (int i = 0; i < QweakSimG4_Tree->GetEntries(); i++) {
    QweakSimG4_Tree->GetEntry(i);
    if(i%10000==1) cout<<"   at event: "<<i<<endl;

    evNr = nrEvts + event->Primary.GetPrimaryEventNumber();
    asymPpM = event->Primary.GetAsymDeno();
    asymPmM = event->Primary.GetAsymNomi();
    
    xi=xI[i];
    yi=yI[i];
    zi=zI[i];
    angXi=aXi[i];
    angYi=aYi[i];
    polTi=polI[i];
    
    interaction.clear();
    trackID.clear();
    
    for (int hit = 0; hit < event->Cerenkov.Detector.GetDetectorNbOfHits(); hit++) {	
      if(event->Cerenkov.Detector.GetDetectorID()[hit]!=3) continue;
      int pTypeHit=event->Cerenkov.Detector.GetParticleType()[hit];
      if(abs(pTypeHit)!=11) continue;
      E=event->Cerenkov.Detector.GetTotalEnergy()[hit];
      
      int parentID=event->Cerenkov.Detector.GetParentID()[hit];
      int tID     =event->Cerenkov.Detector.GetParticleID()[hit];
      
      int hasParent(-1),nInt(-1);
      findInt(interaction,trackID,tID,parentID,hasParent,nInt);
      if(nInt!=1 || hasParent==1) continue;

      primary=0;
      if(tID==1 && parentID==0) primary=1;
      double polX=event->Cerenkov.Detector.GetPolarizationX()[hit];
      double polY=event->Cerenkov.Detector.GetPolarizationY()[hit];
      polT=sqrt(pow(polX,2)+pow(polY,2));

      x=event->Cerenkov.Detector.GetDetectorLocalPositionX()[hit];
      y=event->Cerenkov.Detector.GetDetectorLocalPositionY()[hit];
      z=event->Cerenkov.Detector.GetDetectorLocalPositionZ()[hit];
      
      double Gphi   = event->Cerenkov.Detector.GetGlobalPhiAngle()[hit];
      double Gtheta = event->Cerenkov.Detector.GetGlobalThetaAngle()[hit];	
      double _Gtheta=Gtheta/180.*pi;
      double _Gphi=(Gphi+90)/180.*pi; //+90 to account for the offset in the output
      angX = atan2(sin(_Gtheta)*cos(_Gphi),cos(_Gtheta)) * 180.0 / pi;
      angY = atan2(sin(_Gtheta)*sin(_Gphi),cos(_Gtheta)) * 180.0 / pi;

      tout->Fill();
    }//nhit
  }//tree entries
  nrEvts = evNr;
}

void readInitial(string fnm){
  xI.clear();yI.clear();zI.clear();aXi.clear();aYi.clear();polI.clear();

  string posfnm=fnm.substr(0,fnm.rfind("/")+1)+"positionMomentum.in";
  string polfnm=fnm.substr(0,fnm.rfind("/")+1)+"polarization.in";
  ifstream finpos(posfnm.c_str());
  ifstream finpol(polfnm.c_str());
  float tmpx,tmpy,tmpz,tmpax,tmpay,tmppx,tmppy;
  while(finpos>>tmpx>>tmpy>>tmpz>>tmpax>>tmpay){
    finpol>>tmppx>>tmppy;
    xI.push_back(tmpx);
    yI.push_back(tmpy);
    zI.push_back(tmpz);
    aXi.push_back(tmpax);
    aYi.push_back(tmpay);
    polI.push_back(sqrt(pow(tmppx,2)+pow(tmppy,2)));
  }
  finpos.close();
  finpol.close();
}

void findInt(std::vector<int> &inter,std::vector<int> &val, int trackID,int parent, int &hasPar, int &nInt){
  int found=0;
  nInt=-2;
  int findParent=0;
  for(unsigned int i=0;i<val.size();i++){
    if(trackID==val[i]){
      inter[i]++;
      nInt=inter[i];
      found++;
    }
    if(parent==val[i])
      findParent++;
  }
  
  if(findParent) hasPar=1;
  else hasPar=0;
  
  if(!found){
    val.push_back(trackID);
    inter.push_back(1);
    nInt=1;
  }
  if(found>1){
    cout<<"multiple entries for track "<<trackID<<endl;
  }
}
