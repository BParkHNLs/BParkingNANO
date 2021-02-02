#!/bin/bash

#--------------------
# This script launches the nanoAOD production of a given file on slurm
# ${1}:  outdir
# ${2}:  usr
# ${3}:  pl
# ${4}:  tag
# ${5}:  isMC
# ${6}:  isRemote
# ${7}:  filelist
# ${8}:  isResubmission (false if first launch)
#--------------------


workdir="/scratch/"${2}"/"${3}"/job_nj"$SLURM_ARRAY_TASK_ID
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
cp run_nano_hnl_cfg.py $workdir

echo "copying the file list(-s) to workdir"
if [ ${8} == 0 ] ; then
  cp ${7} $workdir/filelist.txt
else # different treatment in case of resubmission
  cp -r ${7}* $workdir
  #rm ${7}*$SLURM_ARRAY_TASK_ID*
fi

# index of the output file
outIdx=$SLURM_ARRAY_TASK_ID
echo "index of the outputfile: "$outIdx

cd $workdir

inputFilename=''
if [ ${8} == 0 ] ; then
  inputFilename=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
else # different treatment in case of resubmission
  inputFilename=$(sed '1!d' *nj$SLURM_ARRAY_TASK_ID.txt)
fi
echo "inputfilename: "$inputFilename

if [ ${5} == 1 ] ; then #isMC

  if [ ${6} == 0 ] ; then  #private MC
    echo "going to run nano step on "$inputFilename 
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFile=$inputFilename outputFile="bparknano.root" isMC=True
    DATE_END=`date +%s`
    echo "finished running nano step"
  else # central MC
    echo "going to run nano step on "$inputFilename
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=True
    DATE_END=`date +%s`
    echo "finished running nano step"
  fi

else #isData

  echo "going to run nano step on "$inputFilename
  DATE_START=`date +%s`
  cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=False
  DATE_END=`date +%s`
  echo "finished running nano step"

fi

echo "copying the file"
if [ ${4} == 0 ] ; then
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_nj$outIdx.root 
else
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$outIdx.root 
fi

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

