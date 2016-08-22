# $Id: GNUmakefile,v 1.5 2006/05/02 01:02:20 grimm Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------


name := QweakSimG4
G4TARGET := $(name)
G4EXLIB := true
G4TMP := tmp

G4DEBUG = 1

G4WORKDIR = `pwd`

.PHONY: all
all: lib bin 

# Standard G4 compilation:
ifndef G4SYSTEM
  G4SYSTEM = Linux-g++
endif
include $(G4INSTALL)/config/binmake.gmk

#-------------
# Adding ROOT Package
#-------------
#
ROOTCFLAGS = `root-config --cflags`
CPPFLAGS  +=  $(ROOTCFLAGS) -fPIC
LDFLAGS   += -Wl,--no-as-needed
LDFLAGS   += `root-config --ldflags --new --libs`

#-------------
# Adding Xerces Package
#-------------
#
CPPFLAGS  += -I$(XERCESCROOT)/include
LDFLAGS   += -L$(XERCESCROOT)/lib

#----------------------------
# Adding OpenCascade Package
#----------------------------
#
#INCFLAGS  += -I$(OPENCASCADE)/include 
#LDFLAGS   += -L$(OPENCASCADE)/lib -lTKSTEP -lTKOffset 
#CPPFLAGS  += -DUSE_NON_CONST -DCSFDB -DNO_CXX_EXCEPTION -DNo_Exception -DNDEBUG -DLIN -DLININTEL -DHAVE_CONFIG_H -DUSE_STL_STREAM -DHAVE_LIMITS

#--------------------------
# Adding hadronic libraries
#--------------------------
#
# define G4LISTS_BASE, if you have your own physics lists area install point
# different from the default one, G4LISTS_BASE should point to the directory
# that contains the subdirectory 'hadronic'.
#
#ifndef G4LISTS_BASE
#  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
#  G4LISTS_BASE = $(G4INSTALL)/physics_lists
#else
#  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
#endif

#-------------------------------------------------------------------------------------
# Before we compile, use rootcint to create the dictionary files for our ROOT classes.
#-------------------------------------------------------------------------------------
# 
$(G4TARGET).cc: src/QweakSimUserMainEvent_Dict.cc

src/QweakSimUserMainEvent_Dict.cc: 
	@echo "Generating dictionary QweakSimUserMainEvent_Dict.cc ..."
#	@rootcint -f src/QweakSimUserMainEvent_Dict.cc -c $(ROOTCFLAGS) -I./include QweakSimUserMainEvent.hh QweakSimUserPrimaryEvent.hh QweakSimUserGEM_MainEvent.hh QweakSimUserGEM_SingleGEMEvent.hh QweakSimUserGEM_WirePlaneEvent.hh QweakSimUserHDC_MainEvent.hh QweakSimUserHDC_SingleHDCEvent.hh QweakSimUserHDC_WirePlaneEvent.hh QweakSimUserVDC_MainEvent.hh QweakSimUserVDC_Config.hh QweakSimUserVDC_DriftCellEvent.hh QweakSimUserVDC_SingleVDCEvent.hh QweakSimUserVDC_WirePlaneEvent.hh QweakSimUserTriggerScintillator_MainEvent.hh QweakSimUserTriggerScintillator_DetectorEvent.hh QweakSimUserTriggerScintillator_PMTEvent.hh QweakSimUserCerenkov_MainEvent.hh QweakSimUserCerenkov_OctantEvent.hh QweakSimUserCerenkov_DetectorEvent.hh QweakSimUserCerenkov_PMTEvent.hh EPEvent.hh QweakSimUserMainEvent_LinkDef.hh

	@rootcint -f src/QweakSimUserMainEvent_Dict.cc -c $(ROOTCFLAGS) -I./include QweakSimUserMainEvent.hh QweakSimUserPrimaryEvent.hh QweakSimUserGEM_MainEvent.hh QweakSimUserGEM_SingleGEMEvent.hh QweakSimUserGEM_WirePlaneEvent.hh QweakSimUserHDC_MainEvent.hh QweakSimUserHDC_SingleHDCEvent.hh QweakSimUserHDC_WirePlaneEvent.hh QweakSimUserVDC_MainEvent.hh QweakSimUserVDC_Config.hh QweakSimUserVDC_DriftCellEvent.hh QweakSimUserVDC_SingleVDCEvent.hh QweakSimUserVDC_WirePlaneEvent.hh QweakSimUserTriggerScintillator_MainEvent.hh QweakSimUserTriggerScintillator_DetectorEvent.hh QweakSimUserTriggerScintillator_PMTEvent.hh QweakSimUserCerenkov_MainEvent.hh QweakSimUserCerenkov_DetectorEvent.hh QweakSimUserCerenkov_PMTEvent.hh QweakSimUserMainEvent_LinkDef.hh

#-----------------------------------------------------------------
# After G4 programs have compiled, create the ROOT shared library.
#-----------------------------------------------------------------
#
RootShared = ./QweakSim.so
RootObjs   = $(G4TMPDIR)/QweakSimUserMainEvent.o $(G4TMPDIR)/QweakSimUserMainEvent_Dict.o

bin: $(RootShared)
$(RootShared): $(RootObjs)
	@echo "Making ROOT shared library..."
	@$(CXX) -shared -O $(RootObjs) -o $@
	@echo "$@ done"

#-----------------------------------------------------------------------------
# For cleaning the ROOT files, add a command to the standard G4 "make clean":
#-----------------------------------------------------------------------------
#
clean::
	@echo "Deleting old dictionaries, shared libraries, output files, documentation and analysis classes..."
	@rm -f *.so src/*Dict* *MakeClass* *~ src/*~ include/*~
	@rm -rf docs
	@rm -rf tmp 
	@rm -rf *# src/*# include/*.#
