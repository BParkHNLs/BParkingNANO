import os
import sys
import glob
from os import path
import ROOT


class TagAndProbeLauncher(object):
  def __init__(self, is_data, out_label='', version_label='', ds=None, tagnano=None, tagflat=None, categorisation=None):
    self.is_data = is_data
    self.out_label = out_label
    self.version_label = version_label
    self.ds = ds
    self.tagnano = tagnano 
    self.tagflat = tagflat
    self.categorisation = categorisation
    self.user = 'anlyon'
    self.tree_name = 'tree'
    self.do_submit_batch = True
    self.do_resubmit = False
    self.do_short = False
    self.max_events = 5e6 # this corresponds to the number of events to process per job

    datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
    for dataset in self.ds:
      if self.is_data and dataset not in datasets:
        raise RuntimeError('Invalid dataset, Please choose among {}'.format(datasets))

    categories = ['pt_eta', 'pt_dxysig', 'pt_eta_dxysig']
    if self.categorisation not in categories:
      raise RuntimeError('Invalid categorisation. Please choose among {}'.format(categories))
      
    if not self.is_data:
      self.ds = ['incl']


  def getSubdir(self, ds):
    if self.is_data and ds == 'A1': subdir = 'ParkingBPH1_Run2018A'
    if self.is_data and ds == 'A2': subdir = 'ParkingBPH2_Run2018A'
    if self.is_data and ds == 'A3': subdir = 'ParkingBPH3_Run2018A'
    if self.is_data and ds == 'A4': subdir = 'ParkingBPH4_Run2018A'
    if self.is_data and ds == 'A5': subdir = 'ParkingBPH5_Run2018A'
    if self.is_data and ds == 'A6': subdir = 'ParkingBPH6_Run2018A'
    if self.is_data and ds == 'B1': subdir = 'ParkingBPH1_Run2018B'
    if self.is_data and ds == 'B2': subdir = 'ParkingBPH2_Run2018B'
    if self.is_data and ds == 'B3': subdir = 'ParkingBPH3_Run2018B'
    if self.is_data and ds == 'B4': subdir = 'ParkingBPH4_Run2018B'
    if self.is_data and ds == 'B5': subdir = 'ParkingBPH5_Run2018B'
    if self.is_data and ds == 'B6': subdir = 'ParkingBPH6_Run2018B'
    if self.is_data and ds == 'C1': subdir = 'ParkingBPH1_Run2018C'
    if self.is_data and ds == 'C2': subdir = 'ParkingBPH2_Run2018C'
    if self.is_data and ds == 'C3': subdir = 'ParkingBPH3_Run2018C'
    if self.is_data and ds == 'C4': subdir = 'ParkingBPH4_Run2018C'
    if self.is_data and ds == 'C5': subdir = 'ParkingBPH5_Run2018C'
    if self.is_data and ds == 'D1': subdir = 'ParkingBPH1_Run2018D'
    if self.is_data and ds == 'D2': subdir = 'ParkingBPH2_Run2018D'
    if self.is_data and ds == 'D3': subdir = 'ParkingBPH3_Run2018D'
    if self.is_data and ds == 'D4': subdir = 'ParkingBPH4_Run2018D'
    if self.is_data and ds == 'D5': subdir = 'ParkingBPH5_Run2018D'

    if not self.is_data: subdir = 'BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen'

    return subdir


  def getNEvents(self, flat_file):
    f = ROOT.TFile.Open(flat_file[:flat_file.rfind('.root')+5], 'READ')
    tree = f.Get(self.tree_name)
    n_events = tree.GetEntries()
    f.Close()
    return n_events


  def writeFileList(self, filename, flat_files):
    if not path.exists('./files'):
      os.system('mkdir -p ./files')

    f = open(filename, 'w+')

    n_events = 0
    for flat_file in flat_files:
      n_events += self.getNEvents(flat_file)
      f.write(flat_file+'\n')
    
    f.close()

    name = filename[:filename.rfind('.txt')]

    # divide the filelist in chunks containing approximately max_events
    if self.is_data and n_events > self.max_events:
      number_chunks = int(round(n_events / self.max_events))
      command_split = 'split --number=l/{nchunks} {fn}.txt {fn}_ --additional-suffix=.txt'.format(nchunks=number_chunks, fn=name)
      os.system(command_split)
      os.system('rm {}'.format(filename))

    filelists = [f for f in glob.glob('{}*.txt'.format(name))]

    print ' -> {}*.txt created'.format(name)

    return filelists


  def process(self):
    print ' ----                          ----'
    print '       Tag and Probe Launcher      ' 
    print ' ----                          ----'

    print '\n-> Getting the files'
    flat_files = []
    for dataset in self.ds:
      # get all the files 
      subdir = self.getSubdir(dataset)
      if self.is_data:
        indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/data/{}/{}'.format(self.user, self.version_label, subdir)
      else:
        indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/mc_central/{}/{}'.format(self.user, self.version_label, subdir)

      if self.is_data:
        if self.tagnano != None:
          flat_files += [f for f in glob.glob('{}/Chunk*/flat/flat_bparknano_{}_{}.root'.format(indirectory, self.tagnano, self.tagflat))]
        else:
          flat_files += [f for f in glob.glob('{}/Chunk*/flat/flat_bparknano_{}.root'.format(indirectory, self.tagflat))]
      else:
        if self.tagnano != None:
          flat_files += [f for f in glob.glob('{}/merged/flat_bparknano_{}_{}.root'.format(indirectory, self.tagnano, self.tagflat))]
        else:
          flat_files += [f for f in glob.glob('{}/merged/flat_bparknano_{}.root'.format(indirectory, self.tagflat))]

      if len(flat_files) == 0:
        raise RuntimeError('No sample "{}/Chunk*/flat/flat_bparknano_{}_{}.root" was found - Please check'.format(indirectory, self.tagnano, self.tagflat))

    # create file list
    print '\n-> Preparing the file lists'
    filename = 'files/filelist_{}.txt'.format(self.out_label)
    filelists = self.writeFileList(filename, flat_files)

    # get the submission ready
    if not path.exists('./logs/{}'.format(self.out_label)):
      os.system('mkdir -p ./logs/{}'.format(self.out_label))

    if not path.exists('/work/anlyon/tag_and_probe/outfiles/{}'.format(self.out_label)):
      os.system('mkdir -p /work/anlyon/tag_and_probe/outfiles/{}'.format(self.out_label))

    # submit for each list of files
    print '\n-> Submitting...'
    for ilist, filelist in enumerate(filelists):

      if not self.is_data:
        out_suffix = 'incl' 
      elif self.is_data:
        # get the suffix as the total number of events
        n_events = 0
        f = open(filelist)
        for flat_file in f.readlines():
          n_events += self.getNEvents(flat_file)
        out_suffix = 'chunk{}_n{}'.format(ilist, n_events)
      else:
        out_suffix = 'incl' 

      if self.do_resubmit:
        rootfile = ROOT.TNetXNGFile.Open('/work/anlyon/tag_and_probe/outfiles/{}/results_{}_{}.root'.format(self.out_label, self.out_label, out_suffix), 'READ')
        if not rootfile.TestBit(ROOT.TFile.kRecovered): 
          print '-> previous job was successful'
          continue

      if self.do_submit_batch:
        submit_command = 'sbatch -p {queue} --account t3 --mem 7000 -o ./logs/{out}/log_{sufx}.txt -e ./logs/{out}/log_{sufx}.txt --job-name=tag_and_probe_{out}_{sufx} submitter.sh {infile} {out} {sufx} {cat}'.format(
            #queue = 'standard' if not self.do_short else 'short --time 01:00:00',
            queue = 'long' if not self.do_short else 'short --time 01:00:00',
            infile = filelist,
            out = self.out_label,
            sufx = out_suffix,
            cat = self.categorisation,
            )
      else:
        submit_command = 'sh submitter.sh {infile} {out} {sufx} {cat}'.format(
            infile = filelist,
            out = self.out_label,
            sufx = out_suffix,
            cat = self.categorisation,
            )

      print '\n' + submit_command
      os.system(submit_command)

    print ' - Done -'




if __name__ == "__main__":
  
  #is_data = True
  #out_label = 'test_D1_tag_fired_anyBParkHLT_ptetadxysig_max5e6'
  #version_label = 'V10_30Dec21'
  ##datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  #datasets = ['D1']
  #categorisation = 'pt_eta_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_anyBParkHLT'

  #is_data = True
  #out_label = 'test_fullA_V06_tag_fired_HLT_Mu9_IP6_pteta_max3e6'
  #version_label = 'V06_tag_and_probe'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
  #categorisation = 'pt_eta'
  #tagnano = 'tag_and_probe_v2'
  #tagflat = 'fired_HLT_Mu9_IP6'

  #is_data = True
  #out_label = 'test_D1_tag_fired_anyBParkHLT_pteta_max3e6_v3'
  #version_label = 'V10_30Dec21'
  ##datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  ##datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
  #datasets = ['D1']
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_anyBParkHLT_v2'

  #is_data = True
  #out_label = 'test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_max5e6_v2'
  #version_label = 'V10_30Dec21'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  ##datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
  ##datasets = ['D1']
  #categorisation = 'pt_eta_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'test_fullA_tag_fired_HLT_Mu9_IP6_pteta_max4e6_v2'
  #version_label = 'V10_30Dec21'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6'

  #is_data = True
  #out_label = 'test_D1_tag_fired_DST_DoubleMu1_pteta_max3e6'
  #version_label = 'V10_30Dec21'
  #datasets = ['D1']
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'fired_DST_DoubleMu1'

  #is_data = False
  #out_label = 'test_mc_tag_fired_anyBParkHLT_pteta'
  #version_label = 'BToJPsiKstar_V10_30Dec21'
  #datasets = []
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_anyBParkHLT'

  #is_data = False
  #out_label = 'test_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig'
  #version_label = 'BToJPsiKstar_V10_30Dec21'
  #datasets = []
  #categorisation = 'pt_eta_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'test_mc_V06_tag_fired_HLT_Mu9_IP6_pteta'
  #version_label = 'BToJPsiKstar_V0'
  #datasets = []
  #categorisation = 'pt_eta'
  #tagnano = None
  #tagflat = 'fired_HLT_Mu9_IP6'

  #is_data = False
  #out_label = 'test_mc_tag_fired_DST_DoubleMu1_pteta'
  #version_label = 'BToJPsiKstar_V10_30Dec21'
  #datasets = []
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'fired_DST_DoubleMu1'

  #is_data = False
  #out_label = 'test_mc_tag_fired_HLT_Mu12_IP6_pteta_v2'
  #version_label = 'BToJPsiKstar_V10_30Dec21'
  #datasets = []
  #categorisation = 'pt_eta'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'test_fullBPark_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2'
  #version_label = 'V10_30Dec21'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  #categorisation = 'pt_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'test_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_max5e6_v2'
  #version_label = 'V10_30Dec21'
  #datasets = ['D1']
  #categorisation = 'pt_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'test_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig'
  #version_label = 'BToJPsiKstar_V10_30Dec21'
  #datasets = []
  #categorisation = 'pt_dxysig'
  #tagnano = '30Dec21'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'test_V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_newcat_max3e6'
  #version_label = 'V12_08Aug22'
  #datasets = ['D1']
  #categorisation = 'pt_eta_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'test_V12_08Aug22_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptetadxysig_bsrdst_newcat'
  #version_label = 'BToJPsiKstar_V12_08Aug22'
  #datasets = []
  #categorisation = 'pt_eta_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6'
  #version_label = 'V12_08Aug22'
  #datasets = ['D1']
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'V12_08Aug22_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs'
  #version_label = 'BToJPsiKstar_V12_08Aug22'
  #datasets = []
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'test_V12_08Aug22_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_looseid'
  #version_label = 'BToJPsiKstar_V12_08Aug22'
  #datasets = []
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_looseid'

  #is_data = False
  #out_label = 'test_V12_08Aug22_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_scale1p12_smear0p03'
  #version_label = 'BToJPsiKstar_V12_08Aug22'
  #datasets = []
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_scale1p12_smear0p03'

  #is_data = True
  #out_label = 'V12_08Aug22_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_newbinning_v2'
  #version_label = 'V12_08Aug22'
  #datasets = ['D1']
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = False
  #out_label = 'V12_08Aug22_mc_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_newbinning_v2'
  #version_label = 'BToJPsiKstar_V12_08Aug22'
  #datasets = []
  #categorisation = 'pt_dxysig'
  #tagnano = '08Aug22'
  #tagflat = 'tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6'

  #is_data = True
  #out_label = 'V13_06Feb23_D1_ptdxysig_max5e6'
  #version_label = 'V13_06Feb23'
  #datasets = ['D1']
  #categorisation = 'pt_dxysig'
  #tagnano = '06Feb23'
  #tagflat = 'partial'

  #is_data = True
  #out_label = 'V13_06Feb23_fullBpark_ptdxysig_max5e6'
  #version_label = 'V13_06Feb23'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  #categorisation = 'pt_dxysig'
  #tagnano = '06Feb23'
  #tagflat = 'partial'

  #is_data = True
  #out_label = 'V13_06Feb23_fullBpark_ptdxysig_max5e6_v2'
  #version_label = 'V13_06Feb23'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
  #categorisation = 'pt_dxysig'
  #tagnano = '06Feb23'
  #tagflat = 'tag_and_probe'

  #is_data = True
  #out_label = 'V13_06Feb23_A1_6_ptdxysig_max5e6'
  #version_label = 'V13_06Feb23'
  #datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6']
  #categorisation = 'pt_dxysig'
  #tagnano = '06Feb23'
  #tagflat = 'tag_and_probe'

  is_data = True
  out_label = 'V13_06Feb23_A1_ptdxysig_max5e6'
  version_label = 'V13_06Feb23'
  datasets = ['A1']
  categorisation = 'pt_dxysig'
  tagnano = '06Feb23'
  tagflat = 'tag_and_probe'

  TagAndProbeLauncher(is_data=is_data, out_label=out_label, version_label=version_label, ds=datasets, tagnano=tagnano, tagflat=tagflat, categorisation=categorisation).process()





