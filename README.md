# QweakG4DD
Fork of the Qweak G4 simulation 

To get repo:
git clone git@github.com:cipriangal/QweakG4DD

To build:
mkdir build; cd build

cmake ../

make

There are quite a few warnings but overall it should compile without errors (including on Macs).

The two G4 classes added to the repo are the multiple scattering models used by G4 to propage particles in the simulation. For now the way to add the transverse asymmetry to G4 is to replace the existing files (with the same name) in source/processes/electromagnetic/standard/ src and include directories with the ones provided here. Recompilation of G4 is needed afterwards.
