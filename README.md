# QweakG4DD
Fork of the Qweak G4 simulation 

To get repo:
  `git clone git@github.com:cipriangal/QweakG4DD`

To build:
  ```
  mkdir build; cd build
  cmake ../
  make
  ```
There are quite a few warnings but overall it should compile without errors (including on Macs).

Note: 
This package comes with 2 modified physics processes that will replace the standard Multiple Scattering physics processes in G4 (no Geant4 recompilation needed). They can be found in src/processes/electromagnetic. They are turned on automatically when you do cmake and link to a Geant4 version larger than 10.0.0. They can be turned off in the main executable by setting the "USE_CUSTOM_NUCLEAR_SCATTERING" to 0 (in https://github.com/cipriangal/QweakG4DD/blob/master/QweakSimG4.cc#L47). This version of the physics processes was forked at Geant4 version 10.1.p01 and any modifications added to the files is containted between "FIXME" comments.

Their purpose is to either modify the phi diretion (see https://github.com/cipriangal/QweakG4DD/blob/master/src/processes/electromagnetic/standard/QweakSimUrbanMscModel.cc#L975)  of particles to introduce a transverse asymmetry and calculate an asymmetry weight (https://github.com/cipriangal/QweakG4DD/blob/master/src/processes/electromagnetic/standard/QweakSimUrbanMscModel.cc#L988) and propogate it through to the output. 

## Analysis chain:
### Pb propogation
To run the simulation do:
  `build/QweakSimG4 macro/runPMTDD.mac`

This macro loads another macro (https://github.com/cipriangal/QweakG4DD/blob/master/myQweakCerenkovOnly.mac) that sets up the positions. It also has a couple of other important flags inside the file (that should be included if creating your own macro):

- `/PrimaryEvent/SetFixedPosMom true`
  When set to "true" the exec will shoot electrons at a fixed position set by "/PrimaryEvent/SetBeamPosition?" above. This flag can be set to false and the exec will look for two local files: polarization.in positionMomentum.in to read in the initial position, direction and polarization for the primary electron for each event -- for example running 100 ev in this mode needs 2 files with 100 lines each in a specific format. 
  To automate the creation of these two static files there is a root macro available (https://github.com/cipriangal/QweakG4DD/blob/master/rootScripts/samplePrimaryDist.C). One should run the macro before running the main simulation if flag is "false". It takes the number of events, polarization as either "1" or "-1" and octant distribution (for example running: root -l rootScripts\(100,1,204\) will generate the two static files for 100 events with vertical polarization from a main octant distribution on the face of oct 4. 
- `/PrimaryEvent/SetPolarization V`
  Indendently of the position and momentum direction one could choose to not use the polarization static file generated in the macro and set a fixed polarization for the primary electron, however it is recommended that if using "false" for the SetFixedPosMom above one should set this flag to "f". If not using one of the recognized flags it will default to longitudinal polarization (the custom physics libraries are not in use).
- `/PhysicsProcesses/settingFlag 0`
  This is a flag that allows to turn off the track modification from the macro. A reason for doing this would be to just use the weighting calculation (the track modification will cause it to be wrong by "double-dipping").
- `/TrackingAction/TrackingFlag 0`
  This flag allows for speed up of the simulation by only following certain particles (explaination in the macro). Leaving the optical photons in is quite expensive computationally so we mostly run the primaries or the primaries+secondaries and then use a lookup table to convert the electrons that make it to the quartz in a later step (see below).
- `/random/setSeeds largeNr1 largeNr2`
  Set the two seeds needed by Geant4 for the random number generation.

It produces two output files: **QwSim_0.root**(main output) and **o_tuple.root**(depricated -- should be removed at some point). 

### Tree strip
To get the relevant hits at the face of quartz (and strip the rest of the information that is normally not needed for the rest of the analysis) one then uses https://github.com/cipriangal/QweakG4DD/blob/master/QweakTreeBarHits.cc as

  `build/QweakTreeBarHits QwSim_0.root <outputSuffix> [optional: 1 if you ran your simulation with a fixed position]`

or

  `build/QweakTreeBarHits file.list <outputSuffix> [optional: 1 if you ran your simulation with a fixed position]`

where file.list contains a list of paths and output root files from multiple simulations (for example when running on a farm and you have X files with the same configuration -- see above for different flags). The last argument needs to 1 if you ran your simulation from a fixed location on the bar (i.e. not sampled a moustache).

This produces a **o_hits_outputSuffix.root** file that contains:
  - position and energy information at the face of the quartz: "x,y,z,E" 
  - position information at the beginning sim at the face of the Pb: "xi,yi,zi" 
  - flag to denote primary electron: "primary"
  - event number counter: "evNr"
  - lab reference frame phi and theta angles of the hit: "phi, theta"
  - angles along the Y direction (radial) and X direction (along the bar) at the face of the quartz: "angX", "angY" (and their initial values at the face of the Pb "angXi", "angYi")
  - transverse polarization at the face of the Pb (polTi) and quartz (polT)
  - weight factors for the asymmetry calculation asymPM, asymPP 

### Asymmetry calculation
**Track modification analysis**

Using the hit map from the previous step one can apply different optical models to the hits and get PE distributions using https://github.com/cipriangal/QweakG4DD/blob/master/QweakAnaPEs.cc as:

  `buils/anaPE --rootfile <path to o_hits.root file> --barmodel <name>`
  
Note: This should be used for track modification simulations.

The detailed output for o_anaPE\* files can be understood by looking through the QweakAnaPEs source code. In broad strokes:
  - histograms with names *hpe_\** have the PE distribution for all the events
  - histograms with names *hpeAvg_\** have the average PE distribution. This can give you the uncertainty on the mean
  - histograms with names *posPE_\** have the total number of PEs for a particular position 
  - histograms with names *angPE_\** have the total number of PEs for a particular angle (this is the deflection angle)
  - histograms with names *phiPE_\** have the total number of PEs for a particular global phi angle
  - histograms with names *\*AvgProc_\** have the number of PEs similar to the above but for a chunck of data (by default 5%)

All of these histograms are for either the Left or Right tubes (bad nomenclature: for MD3 this is really beam Left or Right) and for Primaries (P), Non-primaries (N), All electrons (A).

For example if we want to calculate the asymmetry as a function of position along the bar for the Left tube we take posPE_L_P from the o_anaPE_V.root (the output from the vertical polarization simulation) and the same histogram from the o_anaPE_mV.root (the output from the minus vertical polarization simulation). For each bin we can calculate an asymmetry as (V-mV)/(V+mV).

**Weighted analysis**

Alternatively if a weighted simulation was generated to get the number of PEs one should use https://github.com/cipriangal/QweakG4DD/blob/master/QweakAnaPEsWgt.cc as:

  `build/anaPEsWgt --rootfile <path to o_hits.root file> --barmodel <name>`
  

**Effective models analysis**

We have the option to run effective models on the hit distribution (asign an asymmetry to each hit based on some criteria -- either position or momentum direction or energy --) using https://github.com/cipriangal/QweakG4DD/blob/master/averageModel.cc as:

  `build/avgModel --rootfile <path to o_hits.root file> --barmodel <name>`

Note: this should ONLY be used on hit distributions that do not have track modification since track modification already puts an asymmetry into the hit distribution. This uses an enhancement of the number of PEs proportional to the asymmetry (hence the name average model). An earlier version of the effective models moved the track or angle and then calculated the number of PEs (https://github.com/cipriangal/QweakG4DD/blob/master/effectiveModel.cc) -- this is not actively worked on right now, since it can take particles in and out of acceptance --.

This produces a "o_avgModel_%s_%s_offset_%d.root" output file (where the %s are different options) that contains histograms as a function of position, angle along the bar with either the number of PEs for the L/R tubes or the asymmetry, or the Abias or DD. It also prints out a table of values for the asymmetries (L/R), DD and Abias for easy copy paste into a spreadsheet or presentation. 


