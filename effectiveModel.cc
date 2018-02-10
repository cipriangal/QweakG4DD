#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>

#include "QweakSimUserMainEvent.hh"
#include "QweakSimSystemOfUnits.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

using namespace std;

const int dimension=5;//3 DoF + 2 PE values
vector<double> scanPoints[dimension];
int debugPrint=0;
vector<double> totPE(4,0);
void readPEs();
void getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
		std::vector<double> points[dimension]);
void getPEs(std::vector<double> in[dimension], std::vector<double> pt,
	    double &outL, double &outR);

const int processModel = 3;
//model,L/R,Upper/Lower
//0= 221nm pos md8Config16 - DA moustache
//1= 3e-4  ang md8Config16 - DA moustache
//2= 4e-6  ang md8Config16 - DA moustache
//3= 4e-6  ang md8Config16 - centered moustache
double asymLimits[4][2][2]={
  {{-0.194,-0.186},{0.200,0.208}},
  {{-16.7,-16.0},{13.,13.4}},
  {{-0.223,-0.216},{0.174,0.182}},
  {{-0.223,-0.216},{0.174,0.182}}
  };

float model(float val,int type);


int main(int argc, char** argv)
{

  if( argc < 2 ) {
    cout<<" usage: build/effModel [root file with hits]"<<endl;
    return 1;
  }

  readPEs();

  string files(argv[1]);

  TFile *fin=TFile::Open(files.c_str(),"READ");
  TTree *t=(TTree*)fin->Get("t");
  int evNr;
  int primary;//0 secondary, 1 primary
  float x,y,z,E;
  float angX,angY;
  float polT;
  t->SetBranchAddress("evNr",&evNr);
  t->SetBranchAddress("primary",&primary);
  t->SetBranchAddress("x",&x);
  t->SetBranchAddress("y",&y);
  t->SetBranchAddress("z",&z);
  t->SetBranchAddress("E",&E);
  t->SetBranchAddress("angX",&angX);
  t->SetBranchAddress("angY",&angY);
  t->SetBranchAddress("polT",&polT);
  
  TFile *fout=new TFile("o_effModel.root","RECREATE");  

  string hnm[3]={"A","P","M"};
  string hti[3]={"avg"," + "," - "};
  string lr[2]={"l","r"};
  
  TH1D *hpe[2][3],*posPE[2][3],*angPE[2][3];

  TH1D *asR=new TH1D("asR","asymmetry Right PMT",
		     200,asymLimits[processModel][0][0],asymLimits[processModel][0][1]);
  TH1D *asL=new TH1D("asL","asymmetry Left PMT" ,
		     200,asymLimits[processModel][1][0],asymLimits[processModel][1][1]);
  cout<<" Limits for the asymmetry histograms :"<<endl;
  cout<<asymLimits[processModel][1][0]<<" <Left < "<<asymLimits[processModel][1][1]<<endl;
  cout<<asymLimits[processModel][0][0]<<" <Right< "<<asymLimits[processModel][0][1]<<endl;
  
  for(int i=0;i<2;i++)
    for(int j=0;j<3;j++){
      hpe[i][j] = new TH1D(Form("h%spe%s",lr[i].c_str(),hnm[j].c_str()),Form("%s #PEs %s",lr[i].c_str(),hti[j].c_str()),
			   500,0,500);
      posPE[i][j] = new TH1D(Form("h%spe%s_pos",lr[i].c_str(),hnm[j].c_str()),
			     Form("%s #PEs %s;position [cm]",lr[i].c_str(),hti[j].c_str()),
			     200,-100,100);
      angPE[i][j] = new TH1D(Form("h%spe%s_ang",lr[i].c_str(),hnm[j].c_str()),
			     Form("%s #PEs %s;angle [deg]",lr[i].c_str(),hti[j].c_str()),
			     180,-90,90);
      
    }

  double procPE[2][3];
  double stepSize=0.2;
  double currentStep=stepSize;
  int nev=t->GetEntries();
  double totLpe[3]={0,0,0},totRpe[3]={0,0,0};
  for(int i=0;i<nev;i++){
    t->GetEntry(i);
    if(i%1000000==1) cout<<" at event: "<<i<<" "<<float(i+1)/nev*100<<"%"<<endl;

    if(i==0){
      for(int ii=0;ii<2;ii++)
	for(int jj=0;jj<3;jj++)
	  procPE[ii][jj]=0;
    }else if(float(i+1)/nev*100>currentStep){
      asL->Fill((procPE[0][1]-procPE[0][2])/(procPE[0][1]+procPE[0][2]) * 1e6);
      asR->Fill((procPE[1][1]-procPE[1][2])/(procPE[1][1]+procPE[1][2]) * 1e6);
      cout<<i<<" L "<<(procPE[0][1]-procPE[0][2])/(procPE[0][1]+procPE[0][2]) * 1e6<<endl;
      cout<<i<<" R "<<(procPE[1][1]-procPE[1][2])/(procPE[1][1]+procPE[1][2]) * 1e6<<endl;

      for(int ii=0;ii<2;ii++)
	for(int jj=0;jj<3;jj++){
	  procPE[ii][jj]=0;
	}
      currentStep+=stepSize;
    }
    
    //if(i>3000000) break;
    
    if(E<3) continue;
    if(E>100) E=100;

    if(processModel==3)
      x-=3.335;//center distribution on bar
    double var[3]={x,E,angX};
    double modVar[3];
    for(int ivar=0;ivar<3;ivar++)
      modVar[ivar]=model(var[ivar],ivar);    
    int addMod[3]={0,1,-1};
    double range[3]={90,100,80};
    
    int outOfRange=0;
    for(int ivar=0;ivar<3 && outOfRange==0 ;ivar++)
      for(int iadd=0;iadd<3 && outOfRange==0 ;iadd++)
	if(fabs(var[ivar]+addMod[iadd]*modVar[ivar])>range[ivar])
	  outOfRange++;

    if(outOfRange) continue;

    for(int ipt=0;ipt<3;ipt++){

      std::vector<double> pt(dimension-2,0);//correct point
      std::vector<double> pts[dimension];

      for(int ivar=0;ivar<3;ivar++)
	pt[ivar] = var[ivar] + addMod[ipt] * modVar[ivar];

      if(debugPrint){
	cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl;
	cout<<i<<" "<<x<<" "<<E<<" "<<angX<<endl;
	cout<<" == "<<var[0]<<" "<<var[1]<<" "<<var[2]<<endl;
	cout<<" "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<endl;
	cout<<" "<<modVar[0]<<" "<<modVar[1]<<" "<<modVar[2]<<endl;
      }

      double lpe(-1),rpe(-1);
      getCorners(0,scanPoints[0].size(),0,pt,pts);
      getPEs(pts,pt,lpe,rpe);

      if(lpe==-1 || rpe==-1 ||
	 std::isnan(lpe) || std::isnan(rpe) ||
	 std::isinf(lpe) || std::isinf(rpe)){
	cout<<"Problem with interpolator! "<<endl;
	cout<<" "<<ipt<<" "<<lpe<<" "<<rpe<<endl;
	exit(1);
      }

      hpe[0][ipt]->Fill(lpe);
      hpe[1][ipt]->Fill(rpe);
      posPE[0][ipt]->Fill(pt[0],lpe);
      posPE[1][ipt]->Fill(pt[0],rpe);
      angPE[0][ipt]->Fill(pt[2],lpe);
      angPE[1][ipt]->Fill(pt[2],rpe);
      totLpe[ipt]+=lpe;
      totRpe[ipt]+=rpe;
      procPE[0][ipt]+=lpe;
      procPE[1][ipt]+=rpe;
    }

  }


  fout->cd();
  for(int i=0;i<2;i++)
    for(int j=0;j<3;j++){      
      hpe[i][j]->Write();
      posPE[i][j]->Write();
      angPE[i][j]->Write();
    }
  asL->Write();
  asR->Write();
  
  cout<<"Left: ";
  for(int i=0;i<3;i++){
    cout<<setprecision(12)<<totLpe[i]<<" ";
  }

  cout<<endl<<"Right: ";
  for(int i=0;i<3;i++){
    cout<<setprecision(12)<<totRpe[i]<<" ";
  }
  cout<<endl;
  fout->Close();
  return 0;
}

//future models go here
float model(float val,int type){
  if(processModel==0 && type==0)
    return 2.21e-5;//position alone (this is DM 221nm HC position difference)
  else if(processModel==1 && type==2)
    return 2.97e-4;//DM estimation of 5.2 urad (3e-4 deg)
  else if((processModel==2 || processModel==3) && type==2)
    return 4e-6;//angle to get ~400ppb DD
  return 0;
}


void readPEs(){
  //ifstream fin("input/idealBar_alongDir_acrossAng0_lightPara.txt");
  ifstream fin("input/md8Config16_alongDir_acrossAng0_lightPara.txt");
  //ifstream fin("input/idealBar_alongDir_acrossAng23_lightPara.txt");
  //ifstream fin("input/md8Config16_alongDir_acrossAng23_lightPara.txt");
  if(!fin.is_open()) {
    cout<<" cannot read file for PE parametrization :macros/yl_md3_angle_scan.txt" <<endl;
    exit(2);
  }
  double x1,x2,x3,x4,x5,x6,x7,x8,x9;  
  string data;
  getline(fin,data);

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
