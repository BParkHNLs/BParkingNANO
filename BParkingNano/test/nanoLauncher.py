import sys 
import os
import os.path
from os import path
import glob
import subprocess
import ROOT

from nanoTools import NanoTools
sys.path.append('../data/samples')
from bparkingdata_samples import bpark_samples
from qcdmuenriched_samples import qcd_samples
from signal_samples_Jun22 import signal_samples


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to launch the nanoAOD tool on top of privately produced miniAOD files', add_help=True)
  parser.add_argument('--pl'      , type=str, dest='pl'          , help='label of the sample file'                                                            , default=None)
  parser.add_argument('--ds'      , type=str, dest='ds'          , help='[data-mccentral] name of the dataset, e.g /ParkingBPH4/Run2018B-05May2019-v2/MINIAOD', default=None)
  parser.add_argument('--tagnano' , type=str, dest='tagnano'     , help='[optional] tag to be added on the outputfile name of the nano sample'                , default=None)
  parser.add_argument('--tagflat' , type=str, dest='tagflat'     , help='[optional] tag to be added on the outputfile name of the flat sample'                , default=None)
  parser.add_argument('--maxfiles', type=str, dest='maxfiles'    , help='[optional] maximum number of files to process'                                       , default=None)
  parser.add_argument('--user'    , type=str, dest='user'        , help='[optional-mcprivate] specify username where the miniAOD files are stored'            , default=os.environ["USER"])
  parser.add_argument('--mcprivate'         , dest='mcprivate'   , help='run the BParking nano tool on a private MC sample'              , action='store_true', default=False)
  parser.add_argument('--mccentral'         , dest='mccentral'   , help='run the BParking nano tool on a central MC sample'              , action='store_true', default=False)
  parser.add_argument('--sigcentral'        , dest='sigcentral'  , help='run the BParking nano tool on a central signal sample'          , action='store_true', default=False)
  parser.add_argument('--sigcrab'           , dest='sigcrab'     , help='run the BParking nano tool on a CRAB signal sample'             , action='store_true', default=False)
  parser.add_argument('--data'              , dest='data'        , help='run the BParking nano tool on a data sample'                    , action='store_true', default=False)
  parser.add_argument('--donano'            , dest='donano'      , help='launch the nano tool on top of the minifile'                    , action='store_true', default=False)
  parser.add_argument('--doflat'            , dest='doflat'      , help='launch the ntupliser on top of the nanofile'                    , action='store_true', default=False)
  parser.add_argument('--domergenano'       , dest='domergenano' , help='[optional] merge the nanofile steps'                            , action='store_true', default=False)
  parser.add_argument('--dosignal'          , dest='dosignal'    , help='run the BToMuMuPi process'                                      , action='store_true', default=False)
  parser.add_argument('--docontrol'         , dest='docontrol'   , help='run the BToKMuMu process'                                       , action='store_true', default=False)
  parser.add_argument('--dohnl'             , dest='dohnl'       , help='run the HNLToMuMuPi process'                                    , action='store_true', default=False)
  parser.add_argument('--dotageprobe'       , dest='dotageprobe' , help='run the JpsiToMuMu process (tag and probe study)'               , action='store_true', default=False)
  parser.add_argument('--doquick'           , dest='doquick'     , help='[optional] run the jobs on the quick partition (t/job<1h)'      , action='store_true', default=False)
  parser.add_argument('--dolong'            , dest='dolong'      , help='[optional] run the jobs on the long  partition (t/job<7days)'   , action='store_true', default=False)
  parser.add_argument('--dosplitflat'       , dest='dosplitflat' , help='[optional] run the dumper with one job per nano file'           , action='store_true', default=False)
  parser.add_argument('--docompile'         , dest='docompile'   , help='[optional] compile the full BParkingNano tool'                  , action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):
  if opt.pl==None:
    raise RuntimeError('Please indicate the production label: for --mcprivate, it has to correspond to the label of the miniAODl (eg. V01_n9000000_njt300)')

  if opt.data==True and opt.ds==None:
    raise RuntimeError('You are running on data, please indicate the dataset with --ds <dataset>')

  if opt.domergenano==True and opt.donano==False:
    command = 'python nanoMerger.py --donano --pl {}'.format(opt.pl)
    if opt.data==True: command += ' --ds {} --data'.format(opt.ds)
    elif opt.mccentral==True: command += ' --ds {} --mccentral'.format(opt.ds)
    elif opt.sigcentral==True: command += ' --ds {} --sigcentral'.format(opt.ds)
    elif opt.mcprivate==True: command += ' --mcprivate'
    if opt.tagnano != None: command += ' --tagnano {}'.format(opt.tagnano)
    raise RuntimeError('This tool is not well suited for processing the merging of the nano step only. Use instead: \n {}'.format(command))

  if opt.donano==False and opt.doflat==False:
    raise RuntimeError('Please indicate if you want to run the nano tool (--donano) and/or the ntupliser (--doflat)')

  if opt.dosignal==False and opt.docontrol==False and opt.dohnl==False and opt.dotageprobe==False:
    raise RuntimeError('Please indicate the process you want to run (--dosignal and/or --docontrol and/or --dohnl and/or --dotageprobe)')

  if opt.mcprivate==False and opt.mccentral==False and opt.sigcentral==False and opt.sigcrab==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or --mcprivate or --mccentral or --sigcentral or --sigcrab to the command line')

  if opt.mcprivate + opt.mccentral + opt.sigcentral + opt.sigcrab + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral or --sigcentral or --sigcrab to the command line')


class NanoLauncher(NanoTools):
  def __init__(self, opt):
    self.prodlabel   = vars(opt)['pl']
    self.ds          = vars(opt)['ds']
    self.tagnano     = vars(opt)['tagnano']
    self.tagflat     = vars(opt)['tagflat']
    self.maxfiles    = vars(opt)['maxfiles']
    self.mcprivate   = vars(opt)['mcprivate']
    self.mccentral   = vars(opt)['mccentral']
    self.sigcentral  = vars(opt)['sigcentral']
    self.sigcrab     = vars(opt)['sigcrab']
    self.data        = vars(opt)['data']
    self.user        = vars(opt)["user"]
    self.donano      = vars(opt)["donano"]
    self.doflat      = vars(opt)["doflat"]
    self.domergenano = vars(opt)["domergenano"]
    self.dosignal    = vars(opt)["dosignal"]
    self.docontrol   = vars(opt)["docontrol"]
    self.dohnl       = vars(opt)["dohnl"]
    self.dotageprobe = vars(opt)["dotageprobe"]
    self.doquick     = vars(opt)["doquick"]
    self.dolong      = vars(opt)["dolong"]
    self.dosplitflat = vars(opt)["dosplitflat"]
    self.docompile   = vars(opt)["docompile"]

    if self.data:
      if self.ds not in bpark_samples.keys():
        raise RuntimeError('Please indicate on which period of the BParking dataset you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, bpark_samples.keys()))
      self.dataset = bpark_samples[self.ds]
    elif self.mccentral:
      if self.ds not in qcd_samples.keys():
        raise RuntimeError('Please indicate on which QCD dataset you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, qcd_samples.keys()))
      self.dataset = qcd_samples[self.ds]
    elif self.sigcentral:
      if self.ds not in signal_samples.keys():
        raise RuntimeError('Please indicate on which signal dataset you want to run. Label "{}" not recognised. Choose among {}'.format(self.ds, signal_samples.keys()))
      self.dataset = signal_samples[self.ds]


  def compile(self):
    import subprocess
    cwd = os.getcwd()
    loc = cwd[0:cwd.find('PhysicsTools')]
    subprocess.call("scram b", cwd=loc, shell=True)


  def writeFileList(self, nfiles_perchunk, point=None):
    if not path.exists('./files'):
      os.system('mkdir ./files') 

    if self.mcprivate or self.sigcrab:
      filename = './files/filelist_{}_{}'.format(self.prodlabel, point)
    else:
      ds_label = NanoTools.getDataLabel(self, self.dataset) if self.data else NanoTools.getMCLabel(self, self.dataset)
      filename = './files/filelist_{dsl}_{pl}'.format(dsl=ds_label, pl=self.prodlabel) 

    if len(glob.glob('{}*.txt'.format(filename))) == 0: # do not create file if already exists
      if self.mcprivate or self.sigcrab: # fetch miniAOD files
        myfile = open(filename + '.txt', "w+")
        nanofiles = NanoTools.getLocalMiniAODFiles(self, user=self.user, prodlabel=self.prodlabel, point=point, iscrab = True if self.sigcrab else False)

        for nanofile in nanofiles:
          if NanoTools.checkLocalFile(self, nanofile):
            myfile.write(nanofile + '\n')
          else:
            print '    could not open {} --> skipping'.format(nanofile)
      
        myfile.close()  
      else: # fetch files on DAS
        command = 'dasgoclient --query="file dataset={ds} | grep file.name" > {fn}.txt'.format(ds=self.dataset, fn=filename)
        os.system(command)

      # slurm cannot deal with too large arrays
      # -> submit job arrays of size 750
      if NanoTools.getNFiles(self, filename + '.txt') > nfiles_perchunk:
        command_split = 'split -l {nfiles} {fn}.txt {fn}_ --additional-suffix=.txt'.format(nfiles=nfiles_perchunk, fn=filename)
        os.system(command_split)
        os.system('rm {fn}.txt'.format(fn=filename))
        
    print '    ---> {}*.txt created'.format(filename)

    return filename 


  def writeDumperStarter(self, nfiles, outputdir, filelist, label):
    nanoname = 'bparknano' if self.tagnano == None else 'bparknano_{}'.format(self.tagnano) 

    f = open(filelist)
    lines = f.readlines()

    event_chain = []
    event_chain.append('TChain* c = new TChain("Events");')
    for iFile in range(1, nfiles+1):
      file_step = NanoTools.getStep(self, lines[iFile-1], iscrab=True if self.sigcrab else False) if (self.mcprivate or self.sigcrab) else iFile
      #file_step = iFile
      event_chain.append('  c->Add("{}/{}_nj{}.root");'.format(outputdir, nanoname, file_step))
    if self.dosignal:                    event_chain.append('  c->Process("BToMuMuPiDumper.C+", outFileName);')
    if self.docontrol:                   event_chain.append('  c->Process("BToKMuMuDumper.C+", outFileName);')
    if self.dohnl:                       event_chain.append('  c->Process("HNLToMuPiDumper.C+", outFileName);')
    if self.dotageprobe:                 event_chain.append('  c->Process("TagAndProbeDumper.C+", outFileName);')
    #if self.mccentral and self.dosignal: event_chain.append('  c->Process("BackgroundSources.C+", outFileName);')
    event_chain = '\n'.join(event_chain)

    run_chain = []
    run_chain.append('TChain* c_run = new TChain("Runs");')
    for iFile in range(1, nfiles+1):
      file_step = NanoTools.getStep(self, lines[iFile-1], iscrab=True if self.sigcrab else False) if (self.mcprivate or self.sigcrab) else iFile
      #file_step = iFile
      run_chain.append('  c_run->Add("{}/{}_nj{}.root");'.format(outputdir, nanoname, file_step))
    run_chain.append('  c_run->Process("NanoRunDumper.C+", outFileName);')
    run_chain = '\n'.join(run_chain)

    if self.data:
      addMC = ''
    elif self.mcprivate or self.sigcentral or self.sigcrab:
      addMC = 'outFileName += "_isSignalMC";'
    elif self.mccentral:
      addMC = 'outFileName += "_isMC";'

    content = [
      '#include "TChain.h"',
      '#include <iostream>',
      '#include "TProof.h"\n',
      'void starter(){',
      '  TString outFileName = "flat_bparknano.root";',
      '  {addMC}'.format(addMC = addMC),
      '  {addevt}'.format(addevt = event_chain),
      '  {addrun}'.format(addrun = '' if self.data else run_chain),
      '}',
    ]
    content = '\n'.join(content)
          
    if not self.dosplitflat:
      starter_name = './files/starter_{}.C'.format(label)
      dumper_starter = open(starter_name, 'w+')
      dumper_starter.write(content)
      dumper_starter.close()
    else:
      for iFile in range(1, nfiles+1):
        file_step = NanoTools.getStep(self, lines[iFile-1], iscrab=True if self.sigcrab else False) if (self.mcprivate or self.sigcrab) else iFile

        event_chain = []
        event_chain.append('TChain* c = new TChain("Events");')
        event_chain.append('  c->Add("{}/{}_nj{}.root");'.format(outputdir, nanoname, file_step))
        if self.dosignal:    event_chain.append('  c->Process("BToMuMuPiDumper.C+", outFileName);')
        if self.docontrol:   event_chain.append('  c->Process("BToKMuMuDumper.C+", outFileName);')
        if self.dohnl:       event_chain.append('  c->Process("HNLToMuPiDumper.C+", outFileName);')
        if self.dotageprobe: event_chain.append('  c->Process("TagAndProbeDumper.C+", outFileName);')
        event_chain = '\n'.join(event_chain)

        run_chain = []
        run_chain.append('TChain* c_run = new TChain("Runs");')
        run_chain.append('  c_run->Add("{}/{}_nj{}.root");'.format(outputdir, nanoname, file_step))
        run_chain.append('  c_run->Process("NanoRunDumper.C+", outFileName);')
        run_chain = '\n'.join(run_chain)

        content = [
          '#include "TChain.h"',
          '#include <iostream>',
          '#include "TProof.h"\n',
          'void starter(){',
          '  TString outFileName = "flat_bparknano.root";',
          '  {addMC}'.format(addMC = addMC),
          '  {addevt}'.format(addevt = event_chain),
          '  {addrun}'.format(addrun = '' if (self.data or self.dotageprobe) else run_chain),
          '}',
        ]
        content = '\n'.join(content)

        starter_name = './files/starter_{}_nj{}.C'.format(label, iFile)
        dumper_starter = open(starter_name, 'w+')
        dumper_starter.write(content)
        dumper_starter.close()


  def writeMergerSubmitter(self, label, filetype, point=None):
    # defining the command
    command = 'python nanoMerger.py --dobatch --pl {pl}'.format(
        pl = self.prodlabel,
        )

    if self.tagnano != None: command += ' --tagnano {}'.format(self.tagnano)
    if self.tagflat != None: command += ' --tagflat {}'.format(self.tagflat)
    if filetype == 'nano': command += ' --donano' 
    else: command += ' --doflat'
    
    if self.mcprivate or self.sigcrab: 
      command += ' --mcprivate'
      command += ' --point {}'.format(point)
    if self.mccentral: command += ' --ds {} --mccentral'.format(self.ds)
    if self.sigcentral: command += ' --ds {} --sigcentral'.format(self.ds)
    if self.data: command += ' --ds {} --data'.format(self.ds)

    if filetype == 'flat' and self.dosplitflat: command += ' --dosplitflat'

    # defining the workdir
    dirlabel = label
    workdir = '/scratch/{usr}/mergingstep/{lbl}'.format(usr=os.environ["USER"], lbl=dirlabel)

    # content of the submitter
    content = [
      '#!/bin/bash',
      'workdir="{wrkdir}"',
      'mkdir -p $workdir',
      'cp nanoTools.py $workdir',
      'cp nanoMerger.py $workdir',
      'cp haddnano.py $workdir',
      'cp -r ../data/samples/*py $workdir',
      'cd $workdir',
      '{cm}',
      'cd -',
      'rm -r $workdir',
    ]
    content = '\n'.join(content).format(
        wrkdir = workdir,
        cm     = command,
        )

    # create file
    submitter_merger = open('submitter_merger.sh', 'w+')
    submitter_merger.write(content)
    submitter_merger.close()


  def launchNano(self, nfiles, outputdir, logdir, filelist, label):
    slurm_options = '-p {part} --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} --time={hh}:00:00'.format(
      part = 'standard' if not self.doquick else 'short', 
      ld = logdir,
      pl = label,
      ar = '1-{}'.format(nfiles),
      hh = 5 if not self.doquick else 1,
      )

    command = 'sbatch {slurm_opt} submitter.sh {outdir} {usr} {pl} {tag} {isMC} {rmt} {lst} 0 {dosig} {doctrl} {dohnl} {dotep}'.format(
      slurm_opt = slurm_options,
      pl        = label,
      outdir    = outputdir,
      usr       = os.environ["USER"], 
      tag       = 0 if self.tagnano == None else self.tagnano,
      isMC      = 0 if self.data else 1,
      rmt       = 0 if (self.mcprivate or self.sigcrab) else 1,
      lst       = filelist,
      dosig     = 1 if self.dosignal else 0, 
      doctrl    = 1 if self.docontrol else 0, 
      dohnl     = 1 if self.dohnl else 0, 
      dotep     = 1 if self.dotageprobe else 0, 
      )

    job = subprocess.check_output(command, shell=True)
    print '\n       ---> {}'.format(job)

    return NanoTools.getJobId(self, job)


  def launchDumper(self, nfiles, outputdir, logdir, filelist, label, jobId):
    self.writeDumperStarter(nfiles, outputdir, filelist, label)
    tag = NanoTools.getTag(self, self.tagnano, self.tagflat)

    slurm_options = '-p {part} --account=t3 -o {ld}/{ln} -e {ld}/{ln} --job-name=dumperstep_{pl} --time={hh}:00:00 {dp}'.format(
      part    = 'standard' if not self.doquick and not self.dolong else ('short' if self.doquick else 'long'), 
      ln      = 'dumperstep.log' if not self.dosplitflat else 'dumperstep_nj%a.log',
      ld      = logdir,
      pl      = label if not self.dosplitflat else label+'%a',
      hh      = 10 if not self.doquick and not self.dolong else (1 if self.doquick else '7-00'),
      dp      = ('--dependency=afterany:{}'.format(jobId) if jobId != -99 else '') if not self.dosplitflat else '--array 1-{}'.format(nfiles),
      )

    command = 'sbatch {slurm_opt} submitter_dumper.sh {outdir} {usr} {pl} {tag} {isMC} {dosig} {doctrl} {dohnl} {dotep} {splt}'.format(
      slurm_opt = slurm_options,
      pl      = label,
      outdir  = outputdir,
      usr     = os.environ["USER"], 
      tag     = tag,
      isMC    = 0 if self.data else 1,
      dosig   = 1 if self.dosignal else 0,
      doctrl  = 1 if self.docontrol else 0,
      dohnl   = 1 if self.dohnl else 0,
      dotep   = 1 if self.dotageprobe else 0,
      splt    = 0 if not self.dosplitflat else 1,
      )

    job_dump = subprocess.check_output(command, shell=True)
    if jobId == -99:
      print '\n       ---> {}'.format(job_dump)
    else:
      print '       ---> (dependency)' 
      print '            {}'.format(job_dump)

    return NanoTools.getJobId(self, job_dump)


  def launchMerger(self, logdir, label, jobIds, filetype, point=None):
    self.writeMergerSubmitter(label, filetype, point)

    slurm_options = '-p {part} --account=t3 -o {ld}/merger{ft}step.log -e {ld}/merger{ft}step.log --job-name=mergerstep_{pl} --time={hh}:00:00 --dependency=afterany:{jobid}'.format(
      part    = 'standard' if not self.doquick else 'short', 
      ld    = logdir,
      ft    = filetype,
      pl    = label,
      hh    = 10 if not self.doquick else 1,
      jobid = NanoTools.getJobIdsList(self, jobIds),
      )

    command_merge = 'sbatch {slurm_opt} submitter_merger.sh'.format(
        slurm_opt = slurm_options,
        )
    
    n_job_merge = subprocess.check_output(command_merge, shell=True)
    print '       ---> (dependency)' 
    print '            {}'.format(n_job_merge)

    os.system('rm submitter_merger.sh')


  def launchingModule(self, point=None, ds_label=None):
    # declaring some useful quantities
    # ids of the launched jobs, needed for dependency of the merger
    nano_jobIds = []
    flat_jobIds = []

    # counter of processed files
    nfiles_tot = 0
        
    # slurm cannot deal with too large arrays, so does haddnano (keep it hardcoded)
    maxfiles_perchunk = 500
    
    print '\n  --> Fetching the files'
    filelistname = self.writeFileList(maxfiles_perchunk, point)

    # loop on the files (containing at most 500 samples) 
    for iFile, filelist in enumerate(glob.glob('{}*.txt'.format(filelistname))):
      if NanoTools.getNFiles(self, filelist) == 0:
        print '        WARNING: no files were found with the corresponding production label'
        print '                 Did you set the correct username using --user <username>?'
        continue

      # enforcing max files limit
      if self.maxfiles != None and nfiles_tot >= int(self.maxfiles): continue 

      # number of files to process in this chunk
      if self.maxfiles == None or int(self.maxfiles)-nfiles_tot > maxfiles_perchunk:
        nfiles = NanoTools.getNFiles(self, filelist)
      else:
        nfiles = int(self.maxfiles)-nfiles_tot
      nfiles_tot += nfiles
      
      # merging step (if any) must happen after nano or flat steps are completed
      if self.maxfiles == None:
        merge_cond = (iFile == len(glob.glob('{}*.txt'.format(filelistname)))-1)
      else:
        merge_cond = (nfiles_tot == int(self.maxfiles))
        
      print '\n  --> Creating output directory'
      outputdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/'.format(os.environ["USER"])
      if self.mcprivate or self.sigcrab:
        outputdir += '{}/{}/nanoFiles/Chunk{}_n{}'.format(self.prodlabel, point, iFile, nfiles)
      else:
        dirname = 'data' if self.data else ('mc_central' if self.mccentral else 'signal_central')
        outputdir += '{}/{}/{}/Chunk{}_n{}/'.format(dirname, self.prodlabel, ds_label, iFile, nfiles)
      if not path.exists(outputdir):
        os.system('mkdir -p {}'.format(outputdir))
      
      print '\n  --> Creating log directory'
      label1 = self.prodlabel if (self.mcprivate or self.sigcrab) else ds_label
      label2 = point if (self.mcprivate or self.sigcrab) else self.prodlabel
      tag = NanoTools.getTag(self, self.tagnano, self.tagflat)
      #logdir = './logs/{}/{}/Chunk{}_n{}'.format(label1, label2, iFile, nfiles) if tag == 0 \
      #         else './logs/{}/{}_{}/Chunk{}_n{}'.format(label1, label2, tag, iFile, nfiles)
      logdir = '/work/{}/logs/{}/{}/Chunk{}_n{}'.format(os.environ["USER"], label1, label2, iFile, nfiles) if tag == 0 \
               else '/work/{}/logs/{}/{}_{}/Chunk{}_n{}'.format(os.environ["USER"], label1, label2, tag, iFile, nfiles)
      if not path.exists(logdir):
        os.system('mkdir -p {}'.format(logdir))

      label = '{}_{}_Chunk{}_n{}'.format(label1, label2 if tag==0 else label2+'_'+tag, iFile, NanoTools.getNFiles(self, filelist))

      nano_jobId = -99

      if self.donano:
        nano_jobId = self.launchNano(nfiles, outputdir, logdir, filelist, label)

        if self.domergenano:
          nano_jobIds.append(nano_jobId)

          if merge_cond:
            self.launchMerger(logdir, label, nano_jobIds, 'nano')

      if self.doflat:
        flat_outputdir = outputdir + '/flat'
        if not path.exists(flat_outputdir):
          os.system('mkdir -p {}'.format(flat_outputdir))

        flat_jobId = self.launchDumper(nfiles, outputdir, logdir, filelist, label, nano_jobId)

        # merging of flat files happens automatically
        flat_jobIds.append(flat_jobId)

        if merge_cond:
          if point == None:
            self.launchMerger(logdir, label, flat_jobIds, 'flat')
          else:
            self.launchMerger(logdir, label, flat_jobIds, 'flat', point)


  def process(self):
    if self.docompile:
      print '-> Compiling'
      self.compile()

    print '\n------------'
    print ' Processing NanoLauncher on production {} '.format(self.prodlabel if (self.mcprivate or self.sigcrab) else self.dataset)
    print ' --> on the batch'
    print '------------'
    
    if self.mcprivate or self.sigcrab:
      locationSE = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/'.format(self.user, self.prodlabel)

      print '\n-> Getting the different mass points'
      pointsdir = NanoTools.getPointDirs(self, locationSE)
      points    = [point[point.rfind('/')+1:len(point)] for point in pointsdir]
      
      # looping over the signal points
      for point in points:
        print '\n-> Processing mass/ctau point: {}'.format(point)

        #if point != 'mass3.0_ctau184.256851021': continue

        self.launchingModule(point=point)

    
    elif self.data or self.mccentral or self.sigcentral:
      dataset_label = NanoTools.getDataLabel(self, self.dataset) if self.data else NanoTools.getMCLabel(self, self.dataset)

      self.launchingModule(ds_label=dataset_label)

    print '\n-> Submission completed' 

    

if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoLauncher(opt).process()




