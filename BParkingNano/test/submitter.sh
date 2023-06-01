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
# ${9}:  dosignal
# ${10}: docontrol
# ${11}: dohnl
# ${12}: doTagAndProbe
# ${13}: globalTag
# ${14}: force
#--------------------

workdir="/scratch/"${2}"/"${3}"/job_nj"${SLURM_JOB_ID}"_"${SLURM_ARRAY_TASK_ID}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
cp run_nano_hnl_cfg.py $workdir

echo "copying the file list(-s) to workdir"
if [ ${8} == 0 ] ; then
  cp ${7} $workdir/filelist.txt
else # different treatment in case of resubmission
  cp -r ${7}* $workdir
  # copy filelist to pnfs
  #xrdcp -r ${7}*$SLURM_ARRAY_TASK_ID* root://t3dcachedb.psi.ch:1094/${1}
  #rm ${7}*$SLURM_ARRAY_TASK_ID*
  xrdcp -r ${7}"_nj"$SLURM_ARRAY_TASK_ID".txt" root://t3dcachedb.psi.ch:1094/${1}
  rm ${7}"_nj"$SLURM_ARRAY_TASK_ID".txt"
fi

cd $workdir

inputFilename=''
if [ ${8} == 0 ] ; then
  inputFilename=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
else # different treatment in case of resubmission
  if [ ${14} == 0 ] ; then
    inputFilename=$(sed '1!d' *nj$SLURM_ARRAY_TASK_ID.txt)
  else
    inputFilename=$(sed '1!d' *nj$SLURM_ARRAY_TASK_ID.txt)
    echo "inputFilename= "$inputFilename
    echo "xrdcp root://xrootd-cms.infn.it/$inputFilename miniaod.root"
    xrdcp root://xrootd-cms.infn.it/$inputFilename miniaod.root
    inputFilename="miniaod.root"
  fi
fi
echo "inputfilename: "$inputFilename

# index of the output file
if [ ${5} == 1 ] ; then #isMC
  outSuffix=$SLURM_ARRAY_TASK_ID".root" 
  #if [ ${6} == 0 ] ; then  #private MC
  #  search='step4_'
  #  start=$inputFilename
  #  end=${start##*$search}
  #  outIdx=$((${#start} - ${#end}))
  #  outSuffix=${inputFilename:outIdx}
  #else
  #  outSuffix=$SLURM_ARRAY_TASK_ID".root" 
  #fi
else # isData
  outSuffix=$SLURM_ARRAY_TASK_ID".root" 
fi
echo "suffix of the outputfile: "$outSuffix


if [ ${5} == 1 ] ; then #isMC

  if [ ${6} == 0 ] ; then  #private MC
    echo "going to run nano step on "$inputFilename 
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFile=$inputFilename outputFile="bparknano.root" isMC=True doSignal=${9} doControl=${10} doHNL=${11} doTagAndProbe=${12} globalTag=${13}
    DATE_END=`date +%s`
    echo "finished running nano step"
  else # central MC
    echo "going to run nano step on "$inputFilename
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=True doSignal=${9} doControl=${10} doHNL=${11} doTagAndProbe=${12} globalTag=${13}
    DATE_END=`date +%s`
    echo "finished running nano step"
  fi

else #isData

  echo "going to run nano step on "$inputFilename
  DATE_START=`date +%s`
  cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=False doSignal=${9} doControl=${10} doHNL=${11} doTagAndProbe=${12} globalTag=${13}
  DATE_END=`date +%s`
  echo "finished running nano step"

fi

echo "copying the file"
if [ ${4} == 0 ] ; then
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_nj$outSuffix 
else
  echo "xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$outSuffix"
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$outSuffix 
fi

# copy filelist to pnfs
xrdcp -f filelist.txt root://t3dcachedb.psi.ch:1094/${1}/filelist_${3}.txt

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

