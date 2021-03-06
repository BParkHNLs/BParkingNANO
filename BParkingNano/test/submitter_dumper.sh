#!/bin/bash

#--------------------
# This script launches the ntuplising tool 
# ${1}:  outdir  
# ${2}:  usr 
# ${3}:  pl 
# ${4}:  tag 
# ${5}:  isMC
#--------------------

if [ ${5} == 1 ] ; then #isMC
  tag=${4}
else
  tag="0"
fi

workdir="/scratch/"${2}"/"${3}"/dumperjob_"${SLURM_JOB_ID}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying ntupliser to workdir"
cp nanoTools.py $workdir
cp ../data/json/golden_2018.json $workdir
cp ../plugins/dumper/utils.C $workdir 
cp ../plugins/dumper/NanoDumper.C $workdir 
cp ../plugins/dumper/NanoDumper.h $workdir 
cp ../plugins/dumper/NanoRunDumper.C $workdir 
cp ../plugins/dumper/NanoRunDumper.h $workdir 

cd $workdir

echo "creating the starter with command: python nanoTools.py --writestarter --outdir ${1} --tag $tag --ismc ${5}"
python nanoTools.py --writestarter --outdir ${1} --tag $tag --ismc ${5} 

# for test 
cp starter.C $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

echo "running the ntupliser on top of the nanofile"
DATE_START_DUMP=`date +%s`
root -l -q -b "starter.C+" 
DATE_END_DUMP=`date +%s`

echo "copying the file"
if [ ${4} == 0 ] ; then
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano.root
else
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_${4}.root
fi

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime_dump=$((DATE_END_DUMP-DATE_START_DUMP))
echo "Wallclock running time: $runtime_dump s"
