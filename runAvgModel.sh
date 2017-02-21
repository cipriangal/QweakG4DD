#!/bin/sh

if [ $1 == "md1" ]
then 
    bar="md1config10_23"
    hitmap="o_hits_sampledMD1_stdStep_noTrackMod_V.root"
elif [ $1 == "md1model2" ]
then 
    bar="md1config16_model2_23"
    hitmap="o_hits_sampledMD2_stdStep_noTrackMod_V.root"
elif [ $1 == "md1model2" ]
then 
    bar="md1_model2_lightGuideMod"
    hitmap="o_hits_sampledMD2_stdStep_noTrackMod_V.root"
elif [ $1 == "md2" ]
then 
    bar="md2config5_23"
    hitmap="o_hits_sampledMD2_stdStep_noTrackMod_V.root"
elif [ $1 == "md2model2" ]
then 
    bar="md2config5_model2_23"
    hitmap="o_hits_sampledMD2_stdStep_noTrackMod_V.root"
elif [ $1 == "md3" ]
then 
    bar="md3config4_23"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "md4" ]
then 
    bar="md4config4_23"
    hitmap="o_hits_sampledMD4_stdStep_noTrackMod_V.root"
elif [ $1 == "md5" ]
then 
    bar="md5config4_23"
    hitmap="o_hits_sampledMD5_stdStep_noTrackMod_V.root"
elif [ $1 == "md6" ]
then 
    bar="md6config3_23"
    hitmap="o_hits_sampledMD6_stdStep_noTrackMod_V.root"
elif [ $1 == "md7" ]
then 
    bar="md7config2_23"
    hitmap="o_hits_sampledMD7_stdStep_noTrackMod_V.root"
elif [ $1 == "md8" ]
then 
    bar="md8config16_23"
    hitmap="o_hits_sampledMD8_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23" ]
then 
    bar="ideal23"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_polish" ]
then 
    bar="ideal23_polish"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_bevel" ]
then 
    bar="ideal23_bevel"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_glue" ]
then 
    bar="ideal23_glue"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_thickdiff" ]
then 
    bar="ideal23_thickdiff"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_RBevelEndcapCentralGlueSideOnly" ]
then 
    bar="ideal23_RBevelEndcapCentralGlueSideOnly"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_RBevelEndcapPMTSideOnly" ]
then 
    bar="ideal23_RBevelEndcapPMTSideOnly"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_RBevelLongAxisOnly" ]
then 
    bar="ideal23_RBevelLongAxisOnly"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_RLG2mmThinner" ]
then 
    bar="ideal23_RLG2mmThinner"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
elif [ $1 == "ideal23_RNoBevel" ]
then 
    bar="ideal23_RNoBevel"
    hitmap="o_hits_sampledMD3_stdStep_noTrackMod_V.root"
fi

./avgModel --rootfile hitmap/$hitmap --barmodel $bar --scan1fct hitmap/o_gp_dd_fold_correctedCor_E3to30.root 0 --Ecut 3 30 --suffix 3to30

./avgModel --rootfile hitmap/$hitmap --barmodel $bar --scan1fct hitmap/o_gp_dd_fold_correctedCor_E3to2000.root 0 --Ecut 3 2000 --suffix 3to2000

./avgModel --rootfile hitmap/$hitmap --barmodel $bar --scan1fct hitmap/o_gp_dd_fold_correctedCor_E30to100.root 0 --Ecut 30 100 --suffix 30to100

./avgModel --rootfile hitmap/$hitmap --barmodel $bar --scan1fct hitmap/o_gp_dd_fold_correctedCor_E100to2000.root 0 --Ecut 100 2000 --suffix 100to2000