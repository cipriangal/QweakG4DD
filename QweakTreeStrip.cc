// Author Ciprian
#include <string.h>
#include "stdio.h"
#include "iostream"
#include <iomanip>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "QweakSimUserMainEvent.hh"
#include "math.h"
#include "stdlib.h"

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2) {
    cout<<" Please specify infile with paths to output QweakSimG4 trees "<<endl;
    return 1;
  }
  
  string fin(argv[1]);
  TChain* QweakSimG4_Tree = new TChain("QweakSimG4_Tree");
  ifstream ifile(fin.c_str());
  char _data[500];
  int scale=0;
  while(ifile>>_data){
    scale+=10000;
    QweakSimG4_Tree->Add(_data);
    cout<<_data<<endl;
  }
  
  //set addresses of leafs
  QweakSimUserMainEvent* event = 0;
  QweakSimG4_Tree->SetBranchAddress("QweakSimUserMainEvent",&event);
  
  //hit
  int evNr;
  int hitNr;
  char pTypeHit;
  double x,y,z;
  double px,py,pz;
  double Gphi,Gtheta;
  double Lphi,Ltheta;
  double vRel;
  double E;
  
  //event
  double LnPMThit,RnPMThit;
  double LEtot,REtot;
  double asymE,asymPMT;
  
  TFile *fout=new TFile("o_anaTree.root","RECREATE");
  TTree *th=new TTree("th","Stripped QweakSimG4 tree for hits");
  th->Branch("evNr",&evNr,"evNr/I");
  th->Branch("hitNr",&hitNr,"hitNr/I");
  th->Branch("pType",&pTypeHit,"pType/C");
  th->Branch("x",&x,"x/F");
  th->Branch("y",&y,"y/F");
  th->Branch("z",&z,"z/F");
  th->Branch("px",&x,"px/F");
  th->Branch("py",&y,"py/F");
  th->Branch("pz",&z,"pz/F");
  th->Branch("Gph",&Gphi,"Gph/F");
  th->Branch("Lph",&Lphi,"Lph/F");
  th->Branch("Gth",&Gtheta,"Gth/F");
  th->Branch("Lth",&Ltheta,"Lth/F");
  th->Branch("vRel",&vRel,"vRel/F");
  th->Branch("E",&E,"E/F");
  
  TTree *te=new TTree("te","Event level info for e#pm");
  te->Branch("LnPMThit",&LnPMThit,"LnPMThit/F");
  te->Branch("RnPMThit",&RnPMThit,"RnPMThit/F");
  te->Branch("LEtot",&LEtot,"LEtot/F");
  te->Branch("REtot",&REtot,"REtot/F");
  te->Branch("asymE",&asymE,"asymE/F");
  te->Branch("asymPMT",&asymPMT,"asymPMT/F");
  
  double cutoff=1./1.458;
  double emass=0.510998910;//MeV

  cout<<" total nr ev: "<<QweakSimG4_Tree->GetEntries()<<endl;
  
  for (int i = 0; i < QweakSimG4_Tree->GetEntries(); i++) {
    QweakSimG4_Tree->GetEntry(i);
    if(i%10000==1) cout<<" at event: "<<i<<endl;
    //loop over hits
    LnPMThit=-2;
    RnPMThit=-2;
    asymE=-2;
    asymPMT=-2;
    LEtot=0;
    REtot=0;
    
    for (int hit = 0; hit < event->Cerenkov.Detector.GetDetectorNbOfHits(); hit++) {
      if(event->Cerenkov.Detector.GetDetectorID()[hit]!=3) continue;
      if(abs(event->Cerenkov.Detector.GetParticleType()[hit])==11) pTypeHit='e';
      else if(abs(event->Cerenkov.Detector.GetParticleType()[hit])==22) pTypeHit='g';
      else continue;
      
      evNr  = i;
      hitNr = hit;
      
      x=event->Cerenkov.Detector.GetDetectorGlobalPositionX()[hit];
      y=event->Cerenkov.Detector.GetDetectorGlobalPositionY()[hit];
      z=event->Cerenkov.Detector.GetDetectorGlobalPositionZ()[hit];
      
      px=event->Cerenkov.Detector.GetGlobalMomentumX()[hit];
      py=event->Cerenkov.Detector.GetGlobalMomentumY()[hit];
      pz=event->Cerenkov.Detector.GetGlobalMomentumZ()[hit];
      
      E=event->Cerenkov.Detector.GetTotalEnergy()[hit];
      
      double _p=sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      double vRel=sqrt(pow(_p,2)/(pow(emass,2)+pow(_p,2)));// relativistic velocity /c which should be compared to 1/refractive index

      Gphi   = event->Cerenkov.Detector.GetGlobalPhiAngle()[hit];
      Gtheta = event->Cerenkov.Detector.GetGlobalThetaAngle()[hit];
      
      Lphi   = event->Cerenkov.Detector.GetLocalPhiAngle()[hit];
      Ltheta = event->Cerenkov.Detector.GetLocalThetaAngle()[hit];
      
      if(vRel>cutoff){
	if(x>0) REtot+=E;
	else LEtot+=E;
      }
      
      th->Fill();
    }
    
    LnPMThit=event->Cerenkov.PMT.GetPMTLeftNbOfHits()[3];//seems to be the same as NbOfPEs
    RnPMThit=event->Cerenkov.PMT.GetPMTRightNbOfHits()[3];
    
    if( (LnPMThit+RnPMThit)>0 ) asymPMT=(LnPMThit-RnPMThit)/(LnPMThit+RnPMThit);
    if( (LEtot + REtot)>0 ) asymE=(LEtot - REtot)/(LEtot + REtot);
    
    te->Fill();
  }
  
  th->Write();
  te->Write();
  fout->Close();
  
  return 0;
}