/*------------------------------------------------------------------------*//*!

 \file QweakSimRoot.cc

 \brief ROOT wrapper with QweakSim functionality

*//*-------------------------------------------------------------------------*/

// ROOT headers
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>

// QweakSim headers
#include "QweakSimRint.hh"

int main(int argc, char** argv)
{
  // Start QweakSim-Root command prompt
  QweakSimRint* qweaksimrint = new QweakSimRint("QweakSim ROOT Analyzer", &argc, argv);
  // Setup include path
  gROOT->ProcessLine(".include include");
  qweaksimrint->Run();
  delete qweaksimrint;
}
