#!/bin/bash

# ${1} infilename
# ${2} out_label
# ${3} out_suffix

homedir=$PWD

workdir='/scratch/anlyon/tagandprobe/'${2}'/'${3}
echo $workdir

mkdir -p $workdir
cp tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py $workdir
cp ${1} $workdir/file.txt
cd $workdir

echo "running script"
DATE_START=`date +%s`
cmsRun tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py inputFiles_load=file.txt outLabel=${2} outSuffix=${3} 
DATE_END=`date +%s`

echo "coyping the files"
#cp -r results_${2}_${3}.root $homedir/outfiles/${2}
mv results_${2}_${3}.root /work/anlyon/tag_and_probe/outfiles/${2}/

cd $homedir
rm -r $workdir

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

echo "End"

