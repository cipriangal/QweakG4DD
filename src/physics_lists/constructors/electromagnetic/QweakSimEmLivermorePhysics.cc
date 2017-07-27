#include "physics_lists/constructors/electromagnetic/QweakSimEmLivermorePhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParticleDefinition.hh"

#include "G4eMultipleScattering.hh"

#include "processes/electromagnetic/standard/QweakSimUrbanMscModel.hh"
#include "processes/electromagnetic/standard/QweakSimWentzelVIModel.hh"


void QweakSimEmLivermorePhysics::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // high energy limit for e+- scattering models
  G4double highEnergyLimit = 100*MeV;
#if G4VERSION_NUMBER > 1021
  auto aParticleIterator=GetParticleIterator();
#endif
  // Add standard EM Processes
  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-" || particleName == "e+") {

      G4eMultipleScattering* msc = new G4eMultipleScattering;
      QweakSimUrbanMscModel* msc1 = new QweakSimUrbanMscModel(asymInfo);
      QweakSimWentzelVIModel* msc2 = new QweakSimWentzelVIModel(asymInfo);
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      ph->RegisterProcess(msc, particle);

    }
  }

  // Fall back
  G4EmLivermorePhysics::ConstructProcess();
}
