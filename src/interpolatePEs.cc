#include "interpolatePEs.hh"

using namespace std;

interpolatePEs::interpolatePEs(string bar):
  barModel(bar),verbosity(0)
{
  readScan();
}

void interpolatePEs::readScan(){

  string path;
  if("ideal0" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng0_lightPara.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng0_lightPara.txt";
  }else if("md8config16_0" == barModel) {
      cout << "Using input/md8Config16_alongDir_acrossAng0_lightPara.txt" << endl;
      path = "input/md8Config16_alongDir_acrossAng0_lightPara.txt";
  }else if("ideal23" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_lightPara.txt";
  }else if("md8config16_23" == barModel) {
      cout << "Using input/md8Config16_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md8Config16_alongDir_acrossAng23_lightPara.txt";
  }else if("md6config3_23" == barModel) {
      cout << "Using input/md6Config3_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md6Config3_alongDir_acrossAng23_lightPara.txt";
  }else if("md7config2_23" == barModel) {
      cout << "Using input/md7Config2_alongDir_acrossAng23_lightPara.txt" << endl;
      path = "input/md7Config2_alongDir_acrossAng23_lightPara.txt";
  }else if("ideal23_polish" == barModel) {
      cout << "Using input/idealBar_alongDir_acrossAng23_Polish0977.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_Polish0977.txt";
  }else if("ideal23_bevel" == barModel) {
      cout << "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt" << endl;
      path = "input/idealBar_alongDir_acrossAng23_Bevel1mmRightBar05mmLeftBar.txt";
  }else{
    cout << "Cannot match your barModel to available list: update list or check name"<<endl;
    exit(2);
  }

  for(int i=0;i<dimension;i++)
    scanPoints[i].clear();

  ifstream fin(path);
  if(!fin.is_open()) {
    cout<<" cannot read file for PE parametrization :macros/yl_md3_angle_scan.txt" <<endl;
    exit(2);
  }

  //setLlimits  
  energyLowLimit  = 3;
  energyHighLimit = 100;

  if("ideal23_bevel" == barModel || "ideal23_polish" == barModel ||
     "md6config3_23" == barModel || "md7config2_23" == barModel ||
     "md8config16_23" == barModel || "ideal23" == barModel ) {
    xPosLowLimit  = -100;
    xPosHighLimit =  100;
    xAngLowLimit  = -89;
    xAngHighLimit =  89;
  } else {
    xPosLowLimit  = -90;
    xPosHighLimit =  90;
    xAngLowLimit  = -80;
    xAngHighLimit =  80;
  }
        
  
  double x1,x2,x3,x4,x5,x6,x7,x8,x9;  
  string data;
  getline(fin, data);
  while(fin>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>x9){
    scanPoints[0].push_back(x1);//position
    scanPoints[1].push_back(x2);//energy
    scanPoints[2].push_back(x3);//angle
    scanPoints[3].push_back(x6);//LPEs
    scanPoints[4].push_back(x8);//RPEs
    if(verbosity)
      cout<<x1<<" "<<x2<<" "<<x3<<" "<<x6<<" "<<x8<<endl;
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
     isnan(lpe) || isnan(rpe) ||
     isinf(lpe) || isinf(rpe)){
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
