#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>

#include "interpolatePEs.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

using namespace std;

void printInfo(TH1D *hl,TH1D *hr);

const int nModels = 8;
//0=                                     
//1= cnst*sgn(angX) for abs(angX)=[20,40]
//2= cnst*angX                           
//3= cnst*sgn(angX)*angX^2               
//4= cnst*angX^3                         
//5= distribution weighted with P
//6= distribution weighted with M
//7= [5]/[6] microscopic distribution

//model,R/L,Upper/Lower
const int rangeTst=0;
double asymLimits[nModels][2][2]={
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-9.000, 9.000},{-9.000, 9.000}},
  {{-9.000, 9.000},{-9.000, 9.000}},
  {{-9.000, 9.000},{-9.000, 9.000}}
};

float model(float val,int type);


int main(int argc, char** argv)
{

  // Print help
  if( argc == 1 || (0 == strcmp("--help", argv[1]))) {
    cout << " usage: build/avgModel [options]" << endl
         << " --rootfile <path to rootfile>" << endl
         << " --barmodel ideal0, ideal23, ideal23_polish, ideal23_bevel, "
         << "md4config4_23, md6config3_23, md7config2_23, md8config16_0 or md8config16_23"
         << endl
         << " --distmodel mirror (omit for as is)" << endl;
    return 1;
  }
  
  // Read in command line paramaters
  TString barModel = "md8config16_23";
  TString distModel = "asIs";
  TString rootfile = "";
  Int_t offset = 0;
  for(Int_t i = 1; i < argc; i++) {
    if(0 == strcmp("--barmodel", argv[i])) {
      barModel = argv[i+1];
    }else if(0 == strcmp("--distmodel", argv[i])) {
      distModel = argv[i+1];
    }else if(0 == strcmp("--rootfile", argv[i])) {
      rootfile = argv[i+1];
    }else if(0 == strcmp("--offset", argv[i])) {
      offset = atoi(argv[i+1]);
    }
  }
  interpolatePEs interpolator(barModel.Data());
  //interpolator.verbosity=1;
  
  // Print out command line paramaters
  cout << "bar model:  " << barModel << endl
       << "distribution model:  " << distModel << endl
       << "using rootfile:  " << rootfile << endl
       << "using offset:  " << offset << endl;

  TFile *fin=TFile::Open(rootfile.Data(),"READ");
  TTree *t=(TTree*)fin->Get("t");
  int evNr;
  int primary;//0 secondary, 1 primary
  float x,y,z,E;
  float angX,angY;
  float angXi,angYi;
  float polT;
  double asymPpM(0),asymPmM(0);
  t->SetBranchAddress("evNr",&evNr);
  t->SetBranchAddress("primary",&primary);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("z",&z);
  t->SetBranchAddress("E",&E);
  t->SetBranchAddress("angX",&angX);
  t->SetBranchAddress("angY",&angY);
  t->SetBranchAddress("angXi",&angXi);
  t->SetBranchAddress("angYi",&angYi);
  t->SetBranchAddress("polT",&polT);
  if(t->GetListOfBranches()->FindObject("asymPpM")){
    t->SetBranchAddress("asymPpM",&asymPpM);
    t->SetBranchAddress("asymPmM",&asymPmM);
  }
  
  TFile *fout=new TFile(Form("o_avgModel_%s_%s_offset_%d.root", barModel.Data(),
                             distModel.Data(),offset),"RECREATE");

  string lr[2]={"R","L"};
  TH1D *hpe[2][nModels],*posPE[2][nModels],*angPE[2][nModels];
  TH1D *as[2][nModels];
  
  for(int i=0;i<nModels;i++)
    for(int j=0;j<2;j++){
      as[j][i]=new TH1D(Form("as%s_%d",lr[j].c_str(),i),Form("model %d %s PMT;asymmetry [ppm]",i,lr[j].c_str()),
			400,asymLimits[i][j][0],asymLimits[i][j][1]);      
      hpe[j][i] = new TH1D(Form("pe%s_%d",lr[j].c_str(),i),Form("model %d %s #PEs",i,lr[j].c_str()),
			   500,0,500);
      posPE[j][i] = new TH1D(Form("pe%s_pos_%d",lr[j].c_str(),i),
			     Form("model %d %s #PEs;position [cm]",i,lr[j].c_str()),
			     200,-100,100);
      angPE[j][i] = new TH1D(Form("pe%s_ang_%d",lr[j].c_str(),i),
			     Form("model %d %s #PEs;angle [deg]",i,lr[j].c_str()),
			     240,-120,120);
    }
    
  std::vector<double> avgStepL(nModels,0);
  std::vector<double> avgStepR(nModels,0);
  std::vector<double> lAvgTotPE(nModels,0);
  std::vector<double> rAvgTotPE(nModels,0);
  
  double stepSize=0.2;
  double currentStep=stepSize;

  int nev=t->GetEntries();

  for(int i=0;i<nev;i++){
    t->GetEntry(i);
    if(i%1000000==1) cout<<" at event: "<<i<<" "<<float(i+1)/nev*100<<"%"<<endl;

    if(float(i+1)/nev*100>currentStep){
      for(int imod=1;imod<nModels;imod++){
	if(avgStepR[0]>0 && avgStepL[0]>0){
	  as[0][imod]->Fill( avgStepR[imod]/avgStepR[0]*1e6 );
	  as[1][imod]->Fill( avgStepL[imod]/avgStepL[0]*1e6 );
	  if(rangeTst){
	    cout<<i<<" "<<imod<<" R "<<avgStepR[imod]<<" "<<avgStepR[0]<<" "<<avgStepR[imod]/avgStepR[0]*1e6<<endl;
	    cout<<i<<" "<<imod<<" L "<<avgStepL[imod]<<" "<<avgStepL[0]<<" "<<avgStepL[imod]/avgStepL[0]*1e6<<endl;
	  }
	}
	avgStepL[imod]=0;
	avgStepR[imod]=0;
      }
      avgStepL[0]=0;
      avgStepR[0]=0;
	
      currentStep+=stepSize;
    }
    
    if(i>1000000 && rangeTst) break;

    x += offset;

    float flip(1.);
    if(distModel == "mirror")
      flip=-1.;

    double lpe(-1),rpe(-1);
    if(!interpolator.getPEs(E,flip*x,flip*angX,lpe,rpe)) continue;
    // cout<<" asdf L R "<<lpe<<" "<<rpe<<endl;
    // std::cin.ignore();
    
    //nModels-1 because microscopic is special and done at the end
    for(int imod=0;imod<nModels-1;imod++){

      double asym=model(angX-angXi,imod);
      if(imod>4){
	double P=(asymPpM+asymPmM)/2;
	double M=(asymPpM-asymPmM)/2;
	if(imod==5){
	  asym=P;
	}else if(imod==6)
	  asym=M;

	if(distModel == "mirror" ) asym=-asym;
      }

      avgStepL[imod]+=asym*lpe;
      avgStepR[imod]+=asym*rpe;
      lAvgTotPE[imod]+=asym*lpe;
      rAvgTotPE[imod]+=asym*rpe;
      
      hpe[0][imod]->Fill((1.+asym)*rpe);
      posPE[0][imod]->Fill(x,asym*rpe);
      angPE[0][imod]->Fill(angX-angXi,asym*rpe);
      hpe[1][imod]->Fill((1.+asym)*lpe);
      posPE[1][imod]->Fill(x,asym*lpe);
      angPE[1][imod]->Fill(angX-angXi,asym*lpe);
      
    }        
  }

  int posBin=posPE[0][7]->GetXaxis()->GetNbins();
  int angBin=angPE[0][7]->GetXaxis()->GetNbins();
  for(int j=0;j<2;j++){
    for(int i=1;i<=posBin;i++){
      double a=posPE[j][5]->GetBinContent(i);
      double b=posPE[j][6]->GetBinContent(i);
      if(a+b>0)
	posPE[j][7]->SetBinContent(i, (a-b)/(a+b) );
    }
    for(int i=1;i<=angBin;i++){
      double a=angPE[j][5]->GetBinContent(i);
      double b=angPE[j][6]->GetBinContent(i);
      if(a+b>0)
	angPE[j][7]->SetBinContent(i, (a-b)/(a+b) );
    }
  }
  
  cout<<endl<<"total PE average: A_L A_R DD A_ave A_ave/DD"<<endl;
  for(int imod=1;imod<nModels;imod++)
    cout<<imod<<"\t"<<lAvgTotPE[imod]/lAvgTotPE[0]<<"\t"<<rAvgTotPE[imod]/rAvgTotPE[0]
	<<"\t"<<lAvgTotPE[imod]/lAvgTotPE[0]-rAvgTotPE[imod]/rAvgTotPE[0]
	<<"\t"<<(lAvgTotPE[imod]/lAvgTotPE[0]+rAvgTotPE[imod]/rAvgTotPE[0])/2
	<<"\t"<<
      ((lAvgTotPE[imod]/lAvgTotPE[0]+rAvgTotPE[imod]/rAvgTotPE[0])/2)/
      (lAvgTotPE[imod]/lAvgTotPE[0]-rAvgTotPE[imod]/rAvgTotPE[0])<<endl;
  fout->cd();
  TNamed* tn1;                              
  TNamed* tn2;                              
  TNamed* tn3;                              
  if("ideal0" == barModel) {
      tn1 = new TNamed("bar","ideal bar");
      tn2 = new TNamed("angle","angle 0");
  }
  if("md8config16_0" == barModel) {
      tn1 = new TNamed("bar","md8config16");
      tn2 = new TNamed("angle","angle 0");
  }                                         
  if("md4config4_23" == barModel) {
      tn1 = new TNamed("bar","md4config4");
      tn2 = new TNamed("angle","angle 23");
  }                                         
  if("md6config3_23" == barModel) {
      tn1 = new TNamed("bar","md6config3");
      tn2 = new TNamed("angle","angle 23");
  }                                         
  if("md7config2_23" == barModel) {
      tn1 = new TNamed("bar","md7config2");
      tn2 = new TNamed("angle","angle 23");
  }                                         
  if("ideal23" == barModel) {
      tn1 = new TNamed("bar","ideal bar");
      tn2 = new TNamed("angle","angle 23");
  }
  if("md8config16_23" == barModel) {
      tn1 = new TNamed("bar","md8config16");
      tn2 = new TNamed("angle","angle 23");
  }
  if("ideal23_polish" == barModel) {
      tn1 = new TNamed("bar","ideal bar with polish");
      tn2 = new TNamed("angle","angle 23");
  }
  if("ideal23_bevel" == barModel) {
      tn1 = new TNamed("bar","ideal bar with bevel");
      tn2 = new TNamed("angle","angle 23");
  }
  if("mirror" == distModel) {
	  tn3 = new TNamed("distribution", "mirror");
  }
  else {
	  tn3 = new TNamed("distribution", "as is");
  }
  tn1->Write();                              
  tn2->Write();                              
  tn3->Write();

  cout<<endl<<" average asymmetry histogram results: A_L dA_L A_R dA_R DD dDD A_bias dA_bia A_bias/DD*100"<<endl;
  for(int j=0;j<nModels;j++){      
    for(int i=0;i<2;i++){
      hpe[i][j]->Write();
      posPE[i][j]->Write();
      angPE[i][j]->Write();
      as[i][j]->Write();
    }

    if(j>0){
      cout<<j<<"\t";
      printInfo(as[1][j],as[0][j]);
      if(as[0][j]->GetBinContent(0)>0 || as[0][j]->GetBinContent(as[0][j]->GetXaxis()->GetNbins()+1)>0 ||
	 as[1][j]->GetBinContent(0)>0 || as[1][j]->GetBinContent(as[1][j]->GetXaxis()->GetNbins()+1)>0){
	cout<<"!!!!! overUnder flow: R L: "<<endl;
	cout<<as[0][j]->GetBinContent(0)<<"\t"
	    <<as[0][j]->GetBinContent(as[0][j]->GetXaxis()->GetNbins()+1)<<"\t"
	    <<as[1][j]->GetBinContent(0)<<"\t"
	    <<as[1][j]->GetBinContent(as[1][j]->GetXaxis()->GetNbins()+1)<<endl;
      }
    }
  }

  fout->Close();
  return 0;
}

//models go here
float model(float val,int type){

  if(type==0)
    return 1;  
  else if(type==1 && (abs(val)>=20 && abs(val)<40) )
    return 4e-6*val/abs(val);
  else if(type==2)
    return 4e-8*val;
  else if(type==3)
    return 1.5e-9*val/abs(val)*pow(val,2);
  else if(type==4)
    return 4e-11*pow(val,3);
  else if(type==5)
    return 0; //reserved for microscopic model

  return 0;
}

void printInfo(TH1D *hl,TH1D *hr){
  double al=hl->GetMean();
  double dal=hl->GetMeanError();
  double ar=hr->GetMean();
  double dar=hr->GetMeanError();

  cout<<al<<"\t"<<dal<<"\t"<<ar<<"\t"<<dar<<"\t"
      <<al-ar<<"\t"<<sqrt(dar*dar+dal*dal)<<"\t"
      <<(al+ar)/2<<"\t"<<sqrt(dar*dar+dal*dal)/2<<"\t"
      <<((al+ar)/2)/(al-ar)*100<<endl;
}
