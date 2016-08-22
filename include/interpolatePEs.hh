#ifndef interpolatePEs_hh
#define interpolatePEs_hh 1

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

class interpolatePEs
{
private:
  double xPosLowLimit, xPosHighLimit;
  double xAngLowLimit, xAngHighLimit;
  double energyLowLimit, energyHighLimit;

  static const int dimension=5;//3 DoF + 2 PE values
  std::vector<double> scanPoints[dimension];
  std::string barModel;
  
  void getCorners(int lowerIndex, int upperIndex, int depth, std::vector<double> point,
		  std::vector<double> points[dimension]);

  void getPEs(std::vector<double> in[dimension], std::vector<double> pt,
	      double &outL, double &outR);

public:
  interpolatePEs(std::string bar);
  ~interpolatePEs(){};
  
  int verbosity;
  int getPEs(double E, double x, double angX ,double &outL,double &outR);
  void readScan();

};

#endif
