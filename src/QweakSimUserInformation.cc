#include "QweakSimUserInformation.hh"

#include "G4PhysicsOrderedFreeVector.hh"

QweakSimUserInformation::QweakSimUserInformation()
{
  cerenkovSecondaryParticleCount = 0;
  Initialize();
}

QweakSimUserInformation::~QweakSimUserInformation()
{
}

void QweakSimUserInformation::Print() const
{
  G4cout << G4endl;
  G4cout << "**** Calling QweakSimUserInformation::Print() ****" << G4endl;
  G4cout << "Primary Event #:: " << PrimaryEventNumber << G4endl;
  G4cout << "Reaction type:: " << ReactionType << G4endl;
  G4cout << "Beam pos:: " << fPositionX[0] << "\t" << fPositionY[0] << "\t"<<fPositionZ[0] << G4endl;
  G4cout << "Tgt center z:: " << TargetCenterPositionZ << G4endl;
  G4cout << "Tgt length:: " << TargetLength << G4endl;
  G4cout << "Beam mom:: " << fNormMomentumX[0] << "\t" << fNormMomentumY[0] << G4endl;
  G4cout << "Org ver pos:: " << OriginVertexPositionX <<"\t" << OriginVertexPositionY << "\t" << OriginVertexPositionZ << G4endl;
  G4cout << "Org ver ang:: "<< OriginVertexThetaAngle <<"\t"<< OriginVertexPhiAngle << G4endl;
  G4cout << "Org ver mom dir:: " << OriginVertexMomentumDirectionX <<"\t" << OriginVertexMomentumDirectionY << "\t" << OriginVertexMomentumDirectionZ << G4endl;
  G4cout << "Org ver KE:: " << OriginVertexKineticEnergy << G4endl;
  G4cout << "Org ver TE:: " << OriginVertexTotalEnergy << G4endl;
  G4cout << "**** Exiting QweakSimUserInformation::Print() ****" << G4endl;
}


void QweakSimUserInformation::PrintELoss()
{  // various energy losses at the target
  G4cout << G4endl;
  G4cout << "**** Calling QweakSimUserInformation::PrintEloss() ****" << G4endl;
  G4cout << "Ion Eloss (In,Out,Tot):: "<< dEIonIn <<"\t" << dEIonOut <<"\t" << GetdEIonTot() << G4endl;
  G4cout << "Brem Eloss (In,Out,Tot):: "<< dEBremIn <<"\t" << dEBremOut <<"\t" << GetdEBremTot() << G4endl;
  G4cout << "Msc Eloss (In,Out,Tot):: "<< dEMscIn <<"\t" << dEMscOut <<"\t" << GetdEMscTot() << G4endl;
  G4cout << "Net Eloss (In,Out,Tot):: "<< GetdETotIn() <<"\t" << GetdETotOut() <<"\t" << GetdETot() << G4endl;
  G4cout << "**** Leaving QweakSimUserInformation::PrintEloss() ****" << G4endl;
}

void QweakSimUserInformation::Initialize()
{

  fPositionX.clear();
  fPositionY.clear();
  fPositionZ.clear();
  fNormMomentumX.clear();
  fNormMomentumY.clear();

  //default step size in Pb is the std
  geantStepSizeInPb = -0.1*mm;
  
  //default is center, front of Pb for MD3 normal incidence
  fPositionX.push_back(0.0*mm);
  fPositionY.push_back(3350.0*mm);
  fPositionZ.push_back(5600.0*mm);
  fNormMomentumX.push_back(0.0*mrad);
  fNormMomentumY.push_back(0.0*mrad);
  
  TargetCenterPositionZ = -652.67*cm;

  TargetLength = 34.346*cm;

  TargetUSDummyPositionOffsetZ = -178.956*mm;  // Offset from LH2 center
  TargetDSDummyPositionOffsetZ = 168.68*mm;    // Offset from LH2 center
  // Dummy target Z positions are taken from the June 21, 2012 target survey
  // The average US and average DS positions are used.  

  // Thicknesses taken from 2012 target survey (https://qweak.jlab.org/elog/Target/21)
  //TargetThicknessUSALDummy1    = 0.8812*mm;
  //TargetThicknessUSALDummy2    = 1.7987*mm;
  //TargetThicknessUSALDummy4    = 3.6030*mm;
  //TargetThicknessDSALDummy2    = 1.8637*mm;
  //TargetThicknessDSALDummy4    = 3.6828*mm;
  //TargetThicknessDSALDummy8    = 7.1980*mm;
  //TargetThicknessUSCDummy      = 0.9973*mm;
  //TargetThicknessDSCDummy      = 3.1876*mm;

  //G4double densityLH2 = 0.0708 g*cm^-3;
  //G4double densityAL  = 2.80   g*cm^-3;  // From 2012 target survey
  //G4double densityUSC = 1.70   g*cm^-3;  //
  //G4double densityDSC = 2.205   g*cm^-3; //

  // Molar masses taken from PDG:  
  // J. Beringer et al. (Particle Data Group), Phys. Rev. D86, 010001 (2012).
  //G4double massLH2 = 1.00794    g*mol^-1;
  //G4double massAL  = 26.9815386 g*mol^-1;
  //G4double massC   = 12.0107    g*mol^-1;

  //TargetLuminosityLH2 = CalculateLuminosity(massLH2, densityLH2, TargetLength);
  //TargetLuminosityUSALWindow = CalculateLuminosity(massAL, densityAL, TargetThicknessUSALWindow);
  //TargetLuminosityDSALWindow;
  //TargetLuminosityUSALDummy1;
  //TargetLuminosityUSALDummy2;
  //TargetLuminosityUSALDummy4;
  //TargetLuminosityDSALDummy2;
  //TargetLuminosityDSALDummy4;
  //TargetLuminosityDSALDummy8;
  //TargetLuminosityUSCDummy;
  //TargetLuminosityDSCDummy;

  PrimaryEventNumber = 0;
  ReactionRegion = 0;
  ReactionType = 0;         // assign a number to which kind of reaction,e.g. 1 = elastic ep,
  PDGcode = 0;              // particle data code/number for the primary particle, e.g. 3=electron
  TrackID = 0;
  GlobalTime = 0.0;
  PrimaryQ2 = 0.0;
  CrossSection = 0.0;
  CrossSectionWeight = 0.0;
  CrossSectionBornTotal          = 0.0;
  CrossSectionBornInelastic      = 0.0;
  CrossSectionBornQE             = 0.0;
  CrossSectionRadTotal           = 0.0;
  CrossSectionRadElastic         = 0.0;
  CrossSectionRadQE              = 0.0;
  CrossSectionRadDIS             = 0.0;
  CrossSectionRadTotalIntOnly    = 0.0;
  CrossSectionRadElasticIntOnly  = 0.0;
  CrossSectionRadQEIntOnly       = 0.0;
  CrossSectionRadDISIntOnly      = 0.0;
  CrossSectionRadElasticPeak     = 0.0;
  Asymmetry = 0.0;
  OriginVertexPositionX = 0.0;
  OriginVertexPositionY = 0.0;
  OriginVertexPositionZ = 0.0;
  OriginVertexThetaAngle = 0.0;
  OriginVertexPhiAngle = 0.0;
  OriginVertexMomentumDirectionX = 0.0;
  OriginVertexMomentumDirectionY = 0.0;
  OriginVertexMomentumDirectionZ = 0.0;
  OriginVertexKineticEnergy = 0.0;
  OriginVertexTotalEnergy = 0.0;

  fEPrime_Min                    = 0.059*GeV;
  fEPrime_Max                    = 1.159*GeV;
  fThetaAngle_Min                = 4.0*degree;
  fThetaAngle_Max                = 13.5*degree;
  fPhiAngle_Min                  = -16.0*degree;
  fPhiAngle_Max                  = 16.0*degree;
  fBeamEnergy                    = 1.16*GeV;
  fLuminosity                    = 0.0;
  fPhaseSpace                    = 0.0;

  OriginVertexKinematicQ2        = 0.0;
  OriginVertexKinematicNu        = 0.0;
  OriginVertexKinematicX         = 0.0;
  OriginVertexKinematicW         = 0.0;

  EffectiveKinematicQ2           = 0.0;
  EffectiveKinematicNu           = 0.0;
  EffectiveKinematicX            = 0.0;
  EffectiveKinematicW            = 0.0;

  // various energy losses at the target
  ClearELoss(); 

  CerEngDep                      = 0.0;
  
  leftPMTHitValid                = 0;
  rightPMTHitValid               = 0;
  StoreStepVolumeName("none");
  SetEdgeEventDetected(0);

  if(cerenkovSecondaryParticleCount){
    delete[] cerenkovSecondaryPartOrig;
    delete[] cerenkovSecondaryPartMom;
    delete[] cerenkovSecondaryPartEng;
    delete[] cerenkovSecondaryPartCharge;
  }
  cerenkovSecondaryParticleCount = 0;
  cerenkovSecondaryElectronCount = 0;
  cerenkovSecondaryPhotonCount = 0;

  cerenkovOpticalPhotonCount = 0;
  CerenkovPhotonEnergy.clear();
  CerenkovPhotonEnergy.resize(0);

  cerenkovSecondaryPositronCount = 0;
  cerenkovSecondaryPartOrig      = NULL;
  cerenkovSecondaryPartMom      = NULL;
  cerenkovSecondaryPartEng      = NULL;
  cerenkovSecondaryPartCharge      = NULL;

  G4ThreeVector tmp(1000,1000,1000);
  cerenkovEventExitPos = tmp;
  
//   G4double D753WKBS20_QE[65][2] = {
//     {200.0*nanometer, 0.68},
//     {210.0*nanometer, 3.55},
//     {220.0*nanometer, 7.40},
//     {230.0*nanometer, 10.4},
//     {240.0*nanometer, 14.6},
//     {250.0*nanometer, 17.8},
//     {260.0*nanometer, 20.6},
//     {270.0*nanometer, 22.6},
//     {280.0*nanometer, 22.4},
//     {290.0*nanometer, 21.8},
//     {300.0*nanometer, 21.1},
//     {310.0*nanometer, 20.5},
//     {320.0*nanometer, 19.7},
//     {330.0*nanometer, 19.2},
//     {340.0*nanometer, 18.4},
//     {350.0*nanometer, 18.0},
//     {360.0*nanometer, 18.2},
//     {370.0*nanometer, 18.8},
//     {380.0*nanometer, 18.3},
//     {390.0*nanometer, 17.6},
//     {400.0*nanometer, 17.7},
//     {410.0*nanometer, 17.6},
//     {420.0*nanometer, 17.5},
//     {430.0*nanometer, 17.1},
//     {440.0*nanometer, 16.7},
//     {450.0*nanometer, 15.8},
//     {460.0*nanometer, 15.0},
//     {470.0*nanometer, 14.4},
//     {480.0*nanometer, 13.7},
//     {490.0*nanometer, 13.1},
//     {500.0*nanometer, 12.4},
//     {510.0*nanometer, 11.7},
//     {520.0*nanometer, 11.0},
//     {530.0*nanometer, 10.4},
//     {540.0*nanometer, 9.77},
//     {550.0*nanometer, 9.15},
//     {560.0*nanometer, 8.53},
//     {570.0*nanometer, 7.95},
//     {580.0*nanometer, 7.39},
//     {590.0*nanometer, 6.87},
//     {600.0*nanometer, 6.38},
//     {610.0*nanometer, 5.90},
//     {620.0*nanometer, 5.45},
//     {630.0*nanometer, 5.07},
//     {640.0*nanometer, 4.71},
//     {650.0*nanometer, 4.39},
//     {660.0*nanometer, 4.10},
//     {670.0*nanometer, 3.79},
//     {680.0*nanometer, 3.51},
//     {690.0*nanometer, 3.25},
//     {700.0*nanometer, 2.98},
//     {710.0*nanometer, 2.68},
//     {720.0*nanometer, 2.40},
//     {730.0*nanometer, 2.13},
//     {740.0*nanometer, 1.88},
//     {750.0*nanometer, 1.65},
//     {760.0*nanometer, 1.47},
//     {770.0*nanometer, 1.30},
//     {780.0*nanometer, 1.13},
//     {790.0*nanometer, 0.96},
//     {800.0*nanometer, 0.80},
//     {810.0*nanometer, 0.65},
//     {820.0*nanometer, 0.48},
//     {830.0*nanometer, 0.33},
//     {840.0*nanometer, 0.18}
//   };

// if quantum efficiency is processed in the Cerenkov detector class,
// set it to 100% here
  G4double D753WKBS20_QE[65][2] = {
    {200.0*nanometer, 100},
    {210.0*nanometer, 100},
    {220.0*nanometer, 100},
    {230.0*nanometer, 100},
    {240.0*nanometer, 100},
    {250.0*nanometer, 100},
    {260.0*nanometer, 100},
    {270.0*nanometer, 100},
    {280.0*nanometer, 100},
    {290.0*nanometer, 100},
    {300.0*nanometer, 100},
    {310.0*nanometer, 100},
    {320.0*nanometer, 100},
    {330.0*nanometer, 100},
    {340.0*nanometer, 100},
    {350.0*nanometer, 100},
    {360.0*nanometer, 100},
    {370.0*nanometer, 100},
    {380.0*nanometer, 100},
    {390.0*nanometer, 100},
    {400.0*nanometer, 100},
    {410.0*nanometer, 100},
    {420.0*nanometer, 100},
    {430.0*nanometer, 100},
    {440.0*nanometer, 100},
    {450.0*nanometer, 100},
    {460.0*nanometer, 100},
    {470.0*nanometer, 100},
    {480.0*nanometer, 100},
    {490.0*nanometer, 100},
    {500.0*nanometer, 100},
    {510.0*nanometer, 100},
    {520.0*nanometer, 100},
    {530.0*nanometer, 100},
    {540.0*nanometer, 100},
    {550.0*nanometer, 100},
    {560.0*nanometer, 100},
    {570.0*nanometer, 100},
    {580.0*nanometer, 100},
    {590.0*nanometer, 100},
    {600.0*nanometer, 100},
    {610.0*nanometer, 100},
    {620.0*nanometer, 100},
    {630.0*nanometer, 100},
    {640.0*nanometer, 100},
    {650.0*nanometer, 100},
    {660.0*nanometer, 100},
    {670.0*nanometer, 100},
    {680.0*nanometer, 100},
    {690.0*nanometer, 100},
    {700.0*nanometer, 100},
    {710.0*nanometer, 100},
    {720.0*nanometer, 100},
    {730.0*nanometer, 100},
    {740.0*nanometer, 100},
    {750.0*nanometer, 100},
    {760.0*nanometer, 100},
    {770.0*nanometer, 100},
    {780.0*nanometer, 100},
    {790.0*nanometer, 100},
    {800.0*nanometer, 100},
    {810.0*nanometer, 100},
    {820.0*nanometer, 100},
    {830.0*nanometer, 100},
    {840.0*nanometer, 100}
  };


  G4double XP4572_QE[15][2] = {
    {200.0*nanometer ,  0.0} ,
    {250.0*nanometer ,  0.0} ,
    {280.0*nanometer ,  0.1} ,
    {290.0*nanometer ,  0.3} ,
    {300.0*nanometer ,  3.0} ,
    {315.0*nanometer , 10.0} ,
    {330.0*nanometer , 20.0} ,
    {350.0*nanometer , 27.5} ,
    {400.0*nanometer , 26.0} ,
    {480.0*nanometer , 20.0} ,
    {540.0*nanometer , 10.0} ,
    {590.0*nanometer ,  3.0} ,
    {615.0*nanometer ,  1.0} ,
    {640.0*nanometer ,  0.3} ,
    {660.0*nanometer ,  0.1}
  };

  PMTQE_XP4572      = new G4PhysicsOrderedFreeVector();
  PMTQED753WKBS20   = new G4PhysicsOrderedFreeVector();
  G4double E_value;

  for (G4int kk=0; kk<65 ; kk++) 
    {    
      if(kk < 15 ){
	E_value= 2*pi*hbarc/( XP4572_QE[kk][0] *nanometer);
	G4cout << "E_value " << kk << " = " << E_value << " QE = " << XP4572_QE[kk][1] << G4endl;
	PMTQE_XP4572->InsertValues(E_value, XP4572_QE[kk][1]/100.);
      }
      
      E_value= 2*pi*hbarc/( D753WKBS20_QE[kk][0] *nanometer);
      PMTQED753WKBS20->InsertValues(E_value, D753WKBS20_QE[kk][1]/100.);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//G4double QweakSimUserInformation::CalculateLuminosity(G4double mass, G4double density, G4double length)
//{
//G4double L = 0.0;                        // Luminosity
//G4double N_b = 6.241*10^12 Hz*microamp;  // # of particles in the beam (from definition of ampere)
//G4double N_A = 6.02214129*10^23 * mol^-1;// Avagadro's number

//L = N_b*length*density*N_A/mass;

//return L;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double QweakSimUserInformation::GetNumberOfPhotoelectrons(G4double eng)
{
  //return PMTQE_XP4572->Value(eng);
  G4bool isOutRange;
  return PMTQE_XP4572->GetValue(eng,isOutRange);
}

G4double QweakSimUserInformation::GetNumberOfPhotoelectronsS20(G4double eng)
{
  //return PMTQED753WKBS20->Value(eng);
  G4bool isOutRange;
  //std::cout<<"Eng="<<eng<<", QE="<<PMTQED753WKBS20->GetValue(eng,isOutRange)<<std::endl;
  return PMTQED753WKBS20->GetValue(eng,isOutRange);
}

void QweakSimUserInformation::StoreCerenkovSecondaryParticleInfo(G4ThreeVector ev,
								 G4ThreeVector em,
								 G4double eng,
								 G4double charge) 
{
  G4int cnt = cerenkovSecondaryParticleCount;
  G4ThreeVector *tmp1 = NULL;
  G4ThreeVector *tmp2 = NULL;
  G4double      *tmp3 = NULL;
  G4double      *tmp4 = NULL;
 
  if(cnt){
    tmp1 = new G4ThreeVector[cnt];
    tmp2 = new G4ThreeVector[cnt];
    tmp3 = new G4double[cnt];
    tmp4 = new G4double[cnt];
  }

  for(G4int i = 0; i < cnt; i++){
    tmp1[i] = cerenkovSecondaryPartOrig[i];
    tmp2[i] = cerenkovSecondaryPartMom[i];
    tmp3[i] = cerenkovSecondaryPartEng[i];
    tmp4[i] = cerenkovSecondaryPartCharge[i];
  }

  if(cnt && cerenkovSecondaryPartOrig) delete[] cerenkovSecondaryPartOrig;
  if(cnt && cerenkovSecondaryPartMom) delete[] cerenkovSecondaryPartMom;
  if(cnt && cerenkovSecondaryPartEng) delete[] cerenkovSecondaryPartEng;
  if(cnt && cerenkovSecondaryPartCharge) delete[] cerenkovSecondaryPartCharge;

  cerenkovSecondaryPartOrig = new G4ThreeVector[cnt+1];
  cerenkovSecondaryPartMom = new G4ThreeVector[cnt+1];
  cerenkovSecondaryPartEng = new G4double[cnt+1];
  cerenkovSecondaryPartCharge = new G4double[cnt+1];

  for(G4int i = 0; i < cnt; i++) {
    cerenkovSecondaryPartOrig[i] = tmp1[i];
    cerenkovSecondaryPartMom[i] = tmp2[i];
    cerenkovSecondaryPartEng[i] = tmp3[i];
    cerenkovSecondaryPartCharge[i] = tmp4[i];
  }

  cerenkovSecondaryPartOrig[cnt] = ev;
  cerenkovSecondaryPartMom[cnt] = em;
  cerenkovSecondaryPartEng[cnt] = eng;
  cerenkovSecondaryPartCharge[cnt] = charge;

  if(cnt){
    delete[] tmp1;
    delete[] tmp2;
    delete[] tmp3;
    delete[] tmp4;
  }

  if(charge == -1)   cerenkovSecondaryElectronCount++;  
  if(charge == 0)   cerenkovSecondaryPhotonCount++;  
  if(charge == 1)   cerenkovSecondaryPositronCount++;  
  cerenkovSecondaryParticleCount++;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ThreeVector QweakSimUserInformation::GetCerenkovSecondaryParticleOrigin(G4int idx) 
{
  G4ThreeVector tmp(1000,1000,1000);
  if(!cerenkovSecondaryParticleCount) return tmp;
  if(idx < 0 || idx >= cerenkovSecondaryParticleCount) return tmp;
  return cerenkovSecondaryPartOrig[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4ThreeVector QweakSimUserInformation::GetCerenkovSecondaryParticleMomentum(G4int idx) 
{
  G4ThreeVector tmp(1000,1000,1000);
  if(!cerenkovSecondaryParticleCount) return tmp;
  if(idx < 0 || idx >= cerenkovSecondaryParticleCount) return tmp;
  return cerenkovSecondaryPartMom[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double QweakSimUserInformation::GetCerenkovSecondaryParticleEnergy(G4int idx) 
{
  if(!cerenkovSecondaryParticleCount) return 0;
  if(idx < 0 || idx >= cerenkovSecondaryParticleCount) return 0;
  return cerenkovSecondaryPartEng[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double QweakSimUserInformation::GetCerenkovSecondaryParticleCharge(G4int idx) 
{
  if(!cerenkovSecondaryParticleCount) return 0;
  if(idx < 0 || idx >= cerenkovSecondaryParticleCount) return 0;
  return cerenkovSecondaryPartCharge[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void QweakSimUserInformation::ResetCerenkovSecondaryParticleInfo()
{
  if(cerenkovSecondaryParticleCount){
    delete[] cerenkovSecondaryPartOrig;
    delete[] cerenkovSecondaryPartMom;
    delete[] cerenkovSecondaryPartEng;
    delete[] cerenkovSecondaryPartCharge;
  }
  cerenkovOpticalPhotonCount = 0;
  CerenkovPhotonEnergy.clear();
  CerenkovPhotonEnergy.resize(0);

  cerenkovSecondaryParticleCount = 0;
  cerenkovSecondaryElectronCount = 0;
  cerenkovSecondaryPhotonCount = 0;
  cerenkovSecondaryPositronCount = 0;
  cerenkovSecondaryPartOrig      = NULL;
  cerenkovSecondaryPartMom      = NULL;
  cerenkovSecondaryPartEng      = NULL;  
  cerenkovSecondaryPartCharge   = NULL;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Read position and momentum values from a file
void QweakSimUserInformation::ReadInitialPositionMomentum(){
  G4cout<<" ==== Reading initial positions and momentums for primary particle ===="<<G4endl;
  fPositionX.clear();
  fPositionY.clear();
  fPositionZ.clear();
  fNormMomentumX.clear();
  fNormMomentumY.clear();
  
  //int i=0;
  std::ifstream fin("positionMomentum.in");
  if(!fin){
    G4cerr<<"In file: "<<__PRETTY_FUNCTION__<<" at line "<<__LINE__<<G4endl
	  <<"\tcould not find file positionMomentum.in."<<G4endl
	  <<"\tyou specified \"/PrimaryEvent/SetFixedPosMom false\" in your macro file."<<G4endl
	  <<"\tQuitting!";
    exit(1);
  }
  G4double tmpx,tmpy,tmpz,tmpPx,tmpPy;
  while(fin>>tmpx>>tmpy>>tmpz>>tmpPx>>tmpPy){
    fPositionX.push_back(tmpx*cm);
    fPositionY.push_back(tmpy*cm);
    fPositionZ.push_back(tmpz*cm);
    fNormMomentumX.push_back(tmpPx*deg);
    fNormMomentumY.push_back(tmpPy*deg);
    // if(i<10)
    //   G4cout<<"  "<<tmpx*cm<<" "<<tmpy*cm<<" "<<tmpz*cm<<" "<<tmpPx*deg<<" "<<tmpPy*deg<<G4endl;
    // i++;
  }
  fin.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Read polarization values from a file
//   these should be generated by the same function as the position/momentum script
//   there should be 100% correlation between the two
void QweakSimUserInformation::ReadInitialPolarization(){
  G4cout<<" ==== Reading initial polarization for primary particle ===="<<G4endl;
  fNormPolarizationX.clear();
  fNormPolarizationY.clear();

  //  int i=0;
  
  std::ifstream fin("polarization.in");
  G4double tmpx(0),tmpy(0);
  while(fin>>tmpx>>tmpy){
    fNormPolarizationX.push_back(tmpx);
    fNormPolarizationY.push_back(tmpy);
    // if(i<10)
    //   G4cout<<"  "<<i<<" "<<tmpx<<" "<<tmpy<<G4endl;
    // i++;
  }
  fin.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
