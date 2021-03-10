# nanoAOD producer customized for BParking analysis 

## Installation

Setup the environment
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

Import the BParking modifications on the TransientTracks, the KinematicVertexFitter, the ElectronRegression and GBRForest
```
git cms-merge-topic -u amlyon:BHNLNano
```

Clone CMSBParking branch for ElectronIndentification
```
git clone --single-branch --branch from-CMSSW_10_2_15_2020Sept15 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

To run on CRAB, do
```
git cms-addpkg RecoEgamma/ElectronIdentification
mkdir -p $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
cp $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
```

Add the BParkingNano tool and build everything

```
git clone git@github.com:BParkHNLs/BParkingNANO.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 8
```

### After first installation

```shell
cd CMSSW_10_2_15/src/
cmsenv 
```


## Nano samples production

```
cd PhysicsTools/BParkingNano/test/
```

### Locally
Modify the inputFiles in run_nano_hnl_cfg.py, make sure to compile and do

```
cmsRun run_nano_hnl_cfg.py 
```

Runs by default on MC and over -1 events. The outputfile will be saved locally in the current working directory.

### On the batch
Runs on a slurm-based engine. 

Since the outputfiles will be saved on the Storage Element, do not forget to activate your proxy

```
voms-proxy-init --voms cms --valid 186:00
```

Then do

```
python nanoLauncher.py <options>
```
Options:

* Indicate whether to run on data or mc (central or private)
  * --mcprivate or   
  * --mccentral or
  *  --data       
* --pl <prodLabel> 
  * with --mcprivate:  must correspond to the production label of the miniAOD sample (e.g V15_full) 
  * with --mccentral/data: any production label of your choice
* --ds <dataset>:  to be used with --data or --mccentral only. Datasets listed in data/samples 
* Indicate on which steps to run
  * --donano: launch the nano step
  * --doflat: launch the ntuplising step
  * --domergenano: launch the merging tool automatically after the nano step. Not recommended as doubles the storage space needed
* --user <user>: with --mcprivate only; username where the miniAOD samples are stored
* --tag <tag>: optional, tag to be appended to the rootfile name 
* --maxfiles <maxfiles>: optional, maximum number of files to process
* --doquick: optional run on slurm quick partition (time/job < 1h)
* --docompile: optional, compiles the BParkingNano tool before launching

Examples of usage:
```
python nanoLauncher.py --pl V15_full --user mratti --donano --doflat --mcprivate
```
```
python nanoLauncher.py --pl V01 --ds /ParkingBPH1/Run2018A-05May2019-v1/MINIAOD --donano --doflat --data
```
```
python nanoLauncher.py --pl V01 --ds /QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM --donano --doflat --mccentral
```

If not done at launching, you can merge a posteriori the different nano steps by doing

```
python nanoMerger.py --pl <prodLabel> --ds <dataset> --tag <tag> --<mcprivate/mccentral/data>
```

Note that the production label and tag have to be consistent with those of the nanoAOD production.



Note:

To make contributions to the central code, see intructions in https://github.com/CMSBParking/BParkingNANO

