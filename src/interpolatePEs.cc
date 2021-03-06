#include "interpolatePEs.hh"
#include "TRandom3.h"

using namespace std;

interpolatePEs::interpolatePEs(string bar,int PEuncert):
  barModel(bar),verbosity(0),peUncert(PEuncert)
{
  readScan();
}
void interpolatePEs::setLightMap(string bar){
  barModel=bar;
  for(int i=0;i<dimension;i++)
    scanPoints[i].clear();
  readScan();
}
  
void interpolatePEs::readScan(){

  string path;
  if("md1config10_23" == barModel) {
    cout << "Using input/md1Config10_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md1Config10_alongDir_acrossAng23_lightPara.txt";
  }else if("md1config16_model2_23" == barModel) {
    cout << "Using input/md1Config16_model2_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md1Config16_model2_alongDir_acrossAng23_lightPara.txt";
  }else if("md1config5_model2_23" == barModel) {
    path = "input/md1Config5_model2_alongDir_acrossAng23.txt";
    cout << "Using: "<< path << endl;
  }else if("md1_model2_lightGuideMod" == barModel) {
    path = "input/md1_model2_LightGuideMod_alongDir_acrossAng23_lightPara.txt";
    cout << "Using: "<< path << endl;
  }else if("md2config5_23" == barModel) {
    cout << "Using input/md2Config5_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md2Config5_alongDir_acrossAng23_lightPara.txt";
  }else if("md2config5_model2_23" == barModel) {
    cout << "Using input/md2Config5_model2_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md2Config5_model2_alongDir_acrossAng23_lightPara.txt";
  }else if("md2config3run1par_model2_23" == barModel) {
    path = "input/md2Config3Run1Parameters_model2_alongDir_acrossAng23_lightPara.txt";
    cout << "Using "<<path<<endl;
  }else if("md2config11_model2_23" == barModel) {
    path = "input/md2Config11_model2_alongDir_acrossAng23.txt";
    cout << "Using "<<path<<endl;
  }else if("md3config4_23" == barModel) {
    cout << "Using input/md3Config4_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md3Config4_alongDir_acrossAng23_lightPara.txt";
  }else if("md4config4_23" == barModel) {
    cout << "Using input/md4Config4_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md4Config4_alongDir_acrossAng23_lightPara.txt";
  }else if("md5config4_23" == barModel) {
    cout << "Using input/md5Config4_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md5Config4_alongDir_acrossAng23_lightPara.txt";
  }else if("md6config3_23" == barModel) {
    cout << "Using input/md6Config3_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md6Config3_alongDir_acrossAng23_lightPara.txt";
  }else if("md7config2_23" == barModel) {
    cout << "Using input/md7Config2_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md7Config2_alongDir_acrossAng23_lightPara.txt";
  }else if("md8config16_0" == barModel) {
    cout << "Using input/md8Config16_alongDir_acrossAng0_lightPara.txt" << endl;
    path = "input/md8Config16_alongDir_acrossAng0_lightPara.txt";
  }else if("md8config16_23" == barModel) {
    cout << "Using input/md8Config16_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/md8Config16_alongDir_acrossAng23_lightPara.txt";
  }else if("md8configMG_23" == barModel) {
    path = "input/md8ConfigMG_alongDir_acrossAng23.txt";
    cout << "Using : " << path << endl;
  }else if("ideal0" == barModel) {
    cout << "Using input/idealBar_alongDir_acrossAng0_lightPara.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng0_lightPara.txt";
  }else if("ideal23" == barModel) {
    cout << "Using input/idealBar_alongDir_acrossAng23_lightPara.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng23_lightPara.txt";
  }else if("ideal23_polish" == barModel) {
    cout << "Using input/idealBar_alongDir_acrossAng23_Polish0977.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng23_Polish0977.txt";
  }else if("ideal23_bevel" == barModel) {
    cout << "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt";
  }else if("ideal23_glue" == barModel) {
    cout << "input/idealBar_alongDir_acrossAng23_GlueFilm040.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng23_GlueFilm040.txt";
  }else if("ideal23_thickdiff" == barModel) {
    cout << "input/idealBar_alongDir_acrossAng23_ThickDiff1.5mm.txt" << endl;
    path = "input/idealBar_alongDir_acrossAng23_ThickDiff1.5mm.txt";
  }else if("ideal23_RBevelEndcapCentralGlueSideOnly" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_RBevelEndcapCentralGlueSideOnly.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_RBevelEndcapPMTSideOnly" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_RBevelEndcapPMTSideOnly.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_RBevelLongAxisOnly" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_RBevelLongAxisOnly.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_RLG2mmThinner" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_RLG2mmThinner.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_RNoBevel" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_RNoBevel.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_GlueFilmR040" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_GlueFilmR040.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_PolishR005Decrease" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_PolishR005Decrease.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23_PolishR010Decrease" == barModel) {
    path = "input/idealBar_alongDir_acrossAng23_PolishR010Decrease.txt";
    cout << "Using: "<< path << endl;
  }else if("ideal23v2"==barModel){
    path = "input/idealBar_alongDir_acrossAng23v2.txt";
    cout << "Using: "<<path << endl;
  }else if("ideal23_1bevelBug"==barModel){
    path = "input/idealBar_alongDir_acrossAng23_1bevelBug.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md1"==barModel){
    path = "input/md1_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md2"==barModel){
    path = "input/md2_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md3"==barModel){
    path = "input/md3_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md4"==barModel){
    path = "input/md4_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md5"==barModel){
    path = "input/md5_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md6"==barModel){
    path = "input/md6_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md7"==barModel){
    path = "input/md7_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else if("tracking_md8"==barModel){
    path = "input/md8_trackingParametrization_dummy3D.txt";
    cout << "Using: "<<path << endl;
  }else{
    cout << "Cannot match your barModel to available list: update list or check name"<<endl;
    exit(2);
  }

  for(int i=0;i<dimension;i++)
    scanPoints[i].clear();

  ifstream fin(path.c_str());
  if(!fin.is_open()) {
    cout << " cannot read file for PE parametrizations"
         << path << endl;
    exit(2);
  }

  //setLlimits  
  energyLowLimit  = 3;
  energyHighLimit = 100;

  if("md8config16_0" == barModel) {
    xPosLowLimit  = -90;
    xPosHighLimit =  90;
    xAngLowLimit  = -80;
    xAngHighLimit =  80;
  } else if("tracking_md1" == barModel || "tracking_md2" == barModel ||
	    "tracking_md3" == barModel || "tracking_md4" == barModel ||
	    "tracking_md5" == barModel || "tracking_md6" == barModel ||
	    "tracking_md7" == barModel || "tracking_md8" == barModel){
    xPosLowLimit  = -95;
    xPosHighLimit =  95;
    xAngLowLimit  = -89;
    xAngHighLimit =  89;    
  } else {
    xPosLowLimit  = -100;
    xPosHighLimit =  100;
    xAngLowLimit  = -89;
    xAngHighLimit =  89;
  }
        
  
  double x1,x2,x3,x4,x5,x6,x7,x8,x9;  
  string data;
  getline(fin, data);
  TRandom3 gs(0);

  while(fin>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>x9){
    double lpe(x6),rpe(x8);
    if(peUncert){
      lpe = gs.Gaus(x6,x7);
      rpe = gs.Gaus(x8,x9);
    }
    
    scanPoints[0].push_back(x1);//position
    scanPoints[1].push_back(x2);//energy
    scanPoints[2].push_back(x3);//angle
    scanPoints[3].push_back(lpe);//LPEs
    scanPoints[4].push_back(rpe);//RPEs
    if(verbosity)
      cout<<x1<<" "<<x2<<" "<<x3<<" "<<lpe<<" "<<rpe<<endl;
  }
  
  fin.close();
}

int interpolatePEs::getPEs(double E, double x, double angX ,double &outL,double &outR){

  if( E > 100 ) E = 100;
  
  if((E > energyHighLimit || E < energyLowLimit) ||
     (x > xPosHighLimit || x < xPosLowLimit) ||
     (angX > xAngHighLimit || angX < xAngLowLimit) )
    {
      if(verbosity)
	cout<<__LINE__<<"\t"<<__PRETTY_FUNCTION__<<endl
	    <<"\tYou provided (E,x,angX):"<<E<<"\t"<<x<<"\t"<<angX<<endl
	    <<"\tLimits are E: "<<energyLowLimit<<"\t"<<energyHighLimit<<endl
	    <<"\tLimits are x: "<<xPosLowLimit<<"\t"<<xPosHighLimit<<endl
	    <<"\tLimits are angX: "<<xAngLowLimit<<"\t"<<xAngHighLimit<<endl;
      return 0;
    }

  std::vector<double> pt(dimension-2,0);
  std::vector<double> pts[dimension];

  if(verbosity)
    cout<<__LINE__<<"\t"<<__PRETTY_FUNCTION__<<endl
	<<"\tTrying to interpolate (E,x,angX):"<<E<<"\t"<<x<<"\t"<<angX<<endl;
  pt[0]=x;
  pt[1]=E;
  pt[2]=angX;

  double lpe(-1),rpe(-1);
  getCorners(0,scanPoints[0].size(),0,pt,pts);
  getPEs(pts,pt,lpe,rpe);

  if(lpe<0 || rpe<0 ||
     std::isnan(lpe) || std::isnan(rpe) ||
     std::isinf(lpe) || std::isinf(rpe)){
    cout<<"Problem with interpolator! "<<endl;
    cout<<" "<<lpe<<" "<<rpe<<endl;
    exit(1);
  }

  if(verbosity)
    cout<<"\tObtained L R: "<<lpe<<" "<<rpe<<endl;
  
  outL=lpe;
  outR=rpe;
  return 1;
}

void interpolatePEs::getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
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
  
  if(verbosity)
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

  if(verbosity){
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
    
    if(verbosity){
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

  if(verbosity)
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl<<"start upper : "<<depth<<" "<<point[depth]<<" "<<*start<<" "<<*(stop-1)<<" "
	<<int(start-begin)<<" "<<int(stop-begin)<<endl;
  
  lI = int( upper_bound(start,stop,point[depth]) - begin );
  if( point[depth] == *(stop-1) )
    hI = upperIndex;
  else{
    valLarger=*(lower_bound(start,stop,point[depth]));
    hI = int( upper_bound(start,stop,valLarger) - begin );
  }
  if(verbosity){
    cout<<endl<<__LINE__<<" "<<__PRETTY_FUNCTION__<<endl
	<<depth<<" "<<lowerIndex<<" "<<upperIndex<<" "<<lI<<" "<<hI<<endl;
    cout<<" "<<valSmaller<<" "<<point[depth]<<" "<<valLarger<<" "
	<<scanPoints[depth][lI]<<" "<<scanPoints[depth][hI-1]<<endl;
  }
  
  if( point[depth]!= *start )
    getCorners(lI,hI,nextDepth,point,points);
}

void interpolatePEs::getPEs(std::vector<double> in[dimension], std::vector<double> pt,
			    double &outL, double &outR){

  std::vector<double> dm[dimension];

  if(verbosity){
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

  if(verbosity){
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
