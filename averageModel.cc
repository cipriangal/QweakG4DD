#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

const int dimension=5;//3 DoF + 2 PE values
vector<double> scanPoints[dimension];
int debugPrint=0;
void readPEs(TString);
void getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
		std::vector<double> points[dimension]);
void getPEs(std::vector<double> in[dimension], std::vector<double> pt,
	    double &outL, double &outR);
void printInfo(TH1D *hl,TH1D *hr);

const int nModels = 6;
//0=                                     
//1= cnst*sgn(angX) for abs(angX)=[20,40]
//2= cnst*angX                           
//3= cnst*sgn(angX)*angX^2               
//4= cnst*angX^3                         
//5= microscopic model calculation

//model,R/L,Upper/Lower
const int rangeTst=0;
double asymLimits[nModels][2][2]={
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
  {{-0.350,-0.100},{ 0.100, 0.350}},
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
         << "md6config3_23, md7config2_23, md8config16_0 or md8config16_23"
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
    }
    if(0 == strcmp("--distmodel", argv[i])) {
      distModel = argv[i+1];
    }
    if(0 == strcmp("--rootfile", argv[i])) {
      rootfile = argv[i+1];
    }
    if(0 == strcmp("--offset", argv[i])) {
      offset = atoi(argv[i+1]);
    }
  }

  // Print out command line paramaters
  cout << "bar model:  " << barModel << endl
       << "distribution model:  " << distModel << endl
       << "using rootfile:  " << rootfile << endl
       << "using offset:  " << offset << endl;

  readPEs(barModel);

  string files(rootfile);

  TFile *fin=TFile::Open(files.c_str(),"READ");
  TTree *t=(TTree*)fin->Get("t");
  int evNr;
  int primary;//0 secondary, 1 primary
  float x,y,z,E;
  float angX,angY;
  float angXi,angYi;
  float polT;
  double calcAsym(0);
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
  if(t->GetListOfBranches()->FindObject("calcAsym"))
    t->SetBranchAddress("calcAsym",&calcAsym);
  
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
	as[0][imod]->Fill( avgStepR[imod]/avgStepR[0]*1e6 );
	as[1][imod]->Fill( avgStepL[imod]/avgStepL[0]*1e6 );
	if(rangeTst){
	  cout<<i<<" "<<imod<<" R "<<avgStepR[imod]<<" "<<avgStepR[0]<<" "<<avgStepR[imod]/avgStepR[0]*1e6<<endl;
	  cout<<i<<" "<<imod<<" L "<<avgStepL[imod]<<" "<<avgStepL[0]<<" "<<avgStepL[imod]/avgStepL[0]*1e6<<endl;
	}
	avgStepL[imod]=0;
	avgStepR[imod]=0;
      }
      avgStepL[0]=0;
      avgStepR[0]=0;
	
      currentStep+=stepSize;
    }
    
    if(i>1000000 && rangeTst) break;


    if("mirror" == distModel) {
        x=-x;
        angX=-angX;
        angXi=-angXi;
    }

    x += offset;

    if("ideal23_bevel" == barModel || "ideal23_polish" == barModel ||
       "md6config3_23" == barModel || "md7config2_23" == barModel ||
       "md8config16_23" == barModel) {
        if(abs(x)>100) continue;
        if(abs(angX)>89) continue;
    }
    else {
        if(abs(x)>90) continue;
        if(abs(angX)>80) continue;
    }
    if(E<3) continue;
    if(E>100) E=100;
        
    std::vector<double> pt(dimension-2,0);//correct point
    std::vector<double> pts[dimension];

    pt[0]=x;
    pt[1]=E;
    pt[2]=angX;
    
    double lpe(-1),rpe(-1);
    getCorners(0,scanPoints[0].size(),0,pt,pts);
    getPEs(pts,pt,lpe,rpe);

    if(lpe<0 || rpe<0 ||
       isnan(lpe) || isnan(rpe) ||
       isinf(lpe) || isinf(rpe)){
      cout<<"Problem with interpolator! "<<endl;
      cout<<" "<<lpe<<" "<<rpe<<endl;
      exit(1);
    }

    for(int imod=0;imod<nModels;imod++){
      double asym=model(angX-angXi,imod);
      if(imod==5){
	asym=calcAsym;
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

void readPEs(TString barModel){
  TString path;
  if("ideal0" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng0_lightPara.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng0_lightPara.txt";
  }
  if("md8config16_0" == barModel) {
      cout << "Using input/md8Config16_alongDir_acrossAng0_lightPara.txt" << endl;
      path = "input/md8Config16_alongDir_acrossAng0_lightPara.txt";
  }
  if("ideal23" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_lightPara.txt";
  }
  if("md8config16_23" == barModel) {
      cout << "Using input/md8Config16_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md8Config16_alongDir_acrossAng23_lightPara.txt";
  }
  if("md6config3_23" == barModel) {
      cout << "Using input/md6Config3_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md6Config3_alongDir_acrossAng23_lightPara.txt";
  }
  if("md7config2_23" == barModel) {
      cout << "Using input/md7Config2_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md7Config2_alongDir_acrossAng23_lightPara.txt";
  }
  if("ideal23_polish" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng23_Polish0977.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_Polish0977.txt";
  }
  if("ideal23_bevel" == barModel) {
      cout << "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt";
  }
  ifstream fin(path);
  if(!fin.is_open()) {
    cout<<" cannot read file for PE parametrization :macros/yl_md3_angle_scan.txt" <<endl;
    exit(2);
  }


  double x1,x2,x3,x4,x5,x6,x7,x8,x9;  
  string data;
  getline(fin, data);

  for(int i=0;i<dimension;i++)
    scanPoints[i].clear();
  
  while(fin>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>x9){
    scanPoints[0].push_back(x1);//position
    scanPoints[1].push_back(x2);//energy
    scanPoints[2].push_back(x3);//angle
    scanPoints[3].push_back(x6);//LPEs
    scanPoints[4].push_back(x8);//RPEs
    if(debugPrint)
      cout<<x1<<" "<<x2<<" "<<x3<<" "<<x6<<" "<<x8<<endl;
  }
  
  fin.close();
}

void getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
		std::vector<double> points[dimension]){

  if(lowerIndex==-1 || upperIndex==-1 || lowerIndex>upperIndex){
    cout<<"Problem with index: "<<lowerIndex<<" "<<upperIndex<<endl;
    exit(1);
  }
  if(lowerIndex==upperIndex) return;
  
  int lI(-1),hI(-1);
  double valSmaller(999),valLarger(999);
  int nextDepth=depth+1;
  
  std::vector<double>::iterator begin=scanPoints[depth].begin();
  std::vector<double>::iterator start=begin+lowerIndex;
  std::vector<double>::iterator stop =begin+upperIndex;
  
  if(debugPrint)
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	<<"start upper : "<<point[depth]<<" "<<*start<<" "<<*(stop-1)<<" "
	<<int(start-begin)<<" "<<int(stop-begin)<<endl;

  if( point[depth]== *start)
    lI=lowerIndex;
  else if( point[depth] == *(stop-1) ){
    lI = int( lower_bound(start,stop,point[depth]) - begin );
  }else{    
    valSmaller = *( lower_bound(start,stop,point[depth]) - 1 );
    lI = int( lower_bound(start,stop,valSmaller) - begin );
  }
  
  hI = int( upper_bound(start,stop,point[depth]) - begin );

  if(debugPrint){
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl<<depth<<" "<<lowerIndex<<" "<<upperIndex<<" "<<lI<<" "<<hI<<endl;
    cout<<" "<<valSmaller<<"<"<<point[depth]<<"<"<<valLarger<<endl
	<<" range: "<<scanPoints[depth][lI]<<" "<<scanPoints[depth][hI-1]<<endl
	<<"  depth lI hI "<<depth<<" "<<lI<<" "<<hI<<endl;
  }

  if(depth==dimension-3){

    for(int i=0;i<dimension;i++) {
      points[i].push_back(scanPoints[i][lI]);
      points[i].push_back(scanPoints[i][hI]);
    }
    
    if(debugPrint){
      cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	  <<"End lower: "<<endl;
      for(int i=0;i<dimension;i++) cout<<scanPoints[i][lI]<<" ";
      cout<<endl;
      for(int i=0;i<dimension;i++) cout<<scanPoints[i][hI]<<" ";
      cout<<endl<<endl;
    }
    return;
  }else{
    getCorners(lI,hI,nextDepth,point,points);
  }

  if( point[depth] == *(stop-1) ) return;

  if(debugPrint)
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl<<"start upper : "<<depth<<" "<<point[depth]<<" "<<*start<<" "<<*(stop-1)<<" "
	<<int(start-begin)<<" "<<int(stop-begin)<<endl;
  
  lI = int( upper_bound(start,stop,point[depth]) - begin );
  if( point[depth] == *(stop-1) )
    hI = upperIndex;
  else{
    valLarger=*(lower_bound(start,stop,point[depth]));
    hI = int( upper_bound(start,stop,valLarger) - begin );
  }
  if(debugPrint){
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	<<depth<<" "<<lowerIndex<<" "<<upperIndex<<" "<<lI<<" "<<hI<<endl;
    cout<<" "<<valSmaller<<" "<<point[depth]<<" "<<valLarger<<" "
	<<scanPoints[depth][lI]<<" "<<scanPoints[depth][hI-1]<<endl;
  }
  
  if( point[depth]!= *start )
    getCorners(lI,hI,nextDepth,point,points);
}

void getPEs(std::vector<double> in[dimension], std::vector<double> pt,
	    double &outL, double &outR){

  std::vector<double> dm[dimension];

  if(debugPrint){
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	<<"start "<<in[0].size()<<endl;
    for(unsigned long i=0;i<in[0].size();i++){
      for(int j=0;j<dimension;j++)
	cout<<in[j][i]<<" ";
      cout<<endl;
    }
  }
  
  if(in[0].size()==1){
    outL=in[dimension-2].front();
    outR=in[dimension-1].front();
    return;
  }

  int depth(-1);
  for(int j=0;j<dimension-2;j++)
    if(in[j][0]!=in[j][1])
      depth=j;
  
  for(int i=0;i<dimension;i++){
    dm[i].resize(in[i].size());
    for(unsigned long j=0;j<in[i].size();j++)
      dm[i][j]=in[i][j];
  }

  for(int i=0;i<dimension;i++)
    in[i].resize(dm[i].size()/2);

  for(unsigned long i=0;i<dm[0].size();i+=2){

    double l1=dm[dimension-2][i];
    double r1=dm[dimension-1][i];
    double l2=dm[dimension-2][i+1];
    double r2=dm[dimension-1][i+1];
    double fl=l1+(l2-l1)*(pt[depth]-dm[depth][i])/(dm[depth][i+1]-dm[depth][i]);
    double fr=r1+(r2-r1)*(pt[depth]-dm[depth][i])/(dm[depth][i+1]-dm[depth][i]);
      
    for(int j=0;j<dimension-2;j++)
      if(j!=depth)
	in[j][i/2]=dm[j][i];
      else
	in[j][i/2]=pt[depth];
    in[dimension-2][i/2]=fl;
    in[dimension-1][i/2]=fr;      
  }

  if(debugPrint){
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	<<depth<<endl;
    for(unsigned long i=0;i<in[0].size();i++){
      for(int j=0;j<dimension;j++)
	cout<<in[j][i]<<" ";
      cout<<endl;
    }
  }
  getPEs(in,pt,outL,outR);
}
