#!/bin/bash

# ${1} infilename
# ${2} out_label
# ${3} out_suffix

homedir=$PWD

workdir='/scratch/anlyon/tagandprobe/'${2}'/'${3}
echo $workdir

mkdir -p $workdir
cp tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py $workdir
cd $workdir

echo "running script"
cmsRun tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py infilename=${1} outLabel=${2} outSuffix=${3} 

echo "coyping the files"
#cp -r results_${2}_${3}.root $homedir/outfiles/${2}
cp -r results_${2}_${3}.root /work/anlyon/tag_and_probe/outfiles/${2}

cd $homedir
rm -r $workdir

echo "End"

