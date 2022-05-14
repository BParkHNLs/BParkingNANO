import os
import sys
import glob
from os import path
import ROOT


class TagAndProbeLauncher(object):
  def __init__(self, is_data, out_label='', version_label='', ds=None, tagnano=None, tagflat=None):
    self.is_data = is_data
    self.out_label = out_label
    self.version_label = version_label
    self.ds = ds
    self.tagnano = tagnano 
    self.tagflat = tagflat
    self.user = 'anlyon'
    self.do_submit_batch = True

    datasets = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C1', 'C2', 'C3', 'C4', 'C5', 'D1', 'D2', 'D3', 'D4', 'D5']
    if self.is_data and self.ds not in datasets:
      raise RuntimeError('Invalid dataset, Please choose among [{}]'.format(datasets))

    if self.is_data and self.ds == 'A1': self.subdir = 'ParkingBPH1_Run2018A'
    if self.is_data and self.ds == 'A2': self.subdir = 'ParkingBPH2_Run2018A'
    if self.is_data and self.ds == 'A3': self.subdir = 'ParkingBPH3_Run2018A'
    if self.is_data and self.ds == 'A4': self.subdir = 'ParkingBPH4_Run2018A'
    if self.is_data and self.ds == 'A5': self.subdir = 'ParkingBPH5_Run2018A'
    if self.is_data and self.ds == 'A6': self.subdir = 'ParkingBPH6_Run2018A'
    if self.is_data and self.ds == 'B1': self.subdir = 'ParkingBPH1_Run2018B'
    if self.is_data and self.ds == 'B2': self.subdir = 'ParkingBPH2_Run2018B'
    if self.is_data and self.ds == 'B3': self.subdir = 'ParkingBPH3_Run2018B'
    if self.is_data and self.ds == 'B4': self.subdir = 'ParkingBPH4_Run2018B'
    if self.is_data and self.ds == 'B5': self.subdir = 'ParkingBPH5_Run2018B'
    if self.is_data and self.ds == 'B6': self.subdir = 'ParkingBPH6_Run2018B'
    if self.is_data and self.ds == 'C1': self.subdir = 'ParkingBPH1_Run2018C'
    if self.is_data and self.ds == 'C2': self.subdir = 'ParkingBPH2_Run2018C'
    if self.is_data and self.ds == 'C3': self.subdir = 'ParkingBPH3_Run2018C'
    if self.is_data and self.ds == 'C4': self.subdir = 'ParkingBPH4_Run2018C'
    if self.is_data and self.ds == 'C5': self.subdir = 'ParkingBPH5_Run2018C'
    if self.is_data and self.ds == 'D1': self.subdir = 'ParkingBPH1_Run2018D'
    if self.is_data and self.ds == 'D2': self.subdir = 'ParkingBPH2_Run2018D'
    if self.is_data and self.ds == 'D3': self.subdir = 'ParkingBPH3_Run2018D'
    if self.is_data and self.ds == 'D4': self.subdir = 'ParkingBPH4_Run2018D'
    if self.is_data and self.ds == 'D5': self.subdir = 'ParkingBPH5_Run2018D'

    if not self.is_data: self.subdir = 'BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen'


  def process(self):
    # get all the files 
    if self.is_data:
      indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/data/{}/{}'.format(self.user, self.version_label, self.subdir)
    else:
      indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/mc_central/{}/{}'.format(self.user, self.version_label, self.subdir)

    #flat_files = [f for f in glob.glob('{}/Chunk*/flat/flat_bparknano_{}_{}_nj*.root'.format(indirectory, self.tagnano, self.tagflat))]
    if self.is_data:
      if self.tagnano != None:
        flat_files = [f for f in glob.glob('{}/Chunk*/flat/flat_bparknano_{}_{}.root'.format(indirectory, self.tagnano, self.tagflat))]
      else:
        flat_files = [f for f in glob.glob('{}/Chunk*/flat/flat_bparknano_{}.root'.format(indirectory, self.tagflat))]
    else:
      if self.tagnano != None:
        flat_files = [f for f in glob.glob('{}/merged/flat_bparknano_{}_{}.root'.format(indirectory, self.tagnano, self.tagflat))]
      else:
        flat_files = [f for f in glob.glob('{}/merged/flat_bparknano_{}.root'.format(indirectory, self.tagflat))]

    for flat_file in flat_files:
      #print flat_file
      if self.is_data:
        chunk = flat_file[flat_file.find('Chunk')+5:flat_file.find('_n', flat_file.find('Chunk')+1)]
        #step = flat_file[flat_file.rfind('nj')+2:flat_file.rfind('root')-1] 
        # suffix corresponds to the total number of events
        f = ROOT.TFile.Open(flat_file, 'READ')
        tree = f.Get('tree')
        n_events = tree.GetEntries()
        out_suffix = '{}_chunk{}_n{}'.format(self.ds, chunk, n_events) 
      else:
        out_suffix = 'incl' #chunk{}_n{}'.format(chunk, n_events) 
      #print out_suffix

      #if out_suffix != 'chunk9_nj99': continue

      if not path.exists('./logs/{}'.format(self.out_label)):
        os.system('mkdir -p ./logs/{}'.format(self.out_label))

      if not path.exists('/work/anlyon/tag_and_probe/outfiles/{}'.format(self.out_label)):
        os.system('mkdir -p /work/anlyon/tag_and_probe/outfiles/{}'.format(self.out_label))

      if self.do_submit_batch:
        submit_command = 'sbatch -p standard --account t3 -o ./logs/{out}/log_{sufx}.txt -e ./logs/{out}/log_{sufx}.txt --job-name=tag_and_probe_{out}_{sufx} submitter.sh {infile} {out} {sufx}'.format(
        #submit_command = 'sbatch -p standard --account t3 -o ./logs/log.txt -e ./logs/log.txt --job-name=tag_and_probe_{out}_{sufx} submitter.sh {infile} {out} {sufx}'.format(
            infile = flat_file,
            out = self.out_label,
            sufx = out_suffix,
            )
      else:
        submit_command = 'sh submitter.sh {infile} {out} {sufx}'.format(
            infile = flat_file,
            out = self.out_label,
            sufx = out_suffix,
            )

      os.system(submit_command)

      #print submit_command





if __name__ == "__main__":
  
  is_data = True
  out_label = 'scale_factors_tag_fired_anyBParkHLT_mc_v2'
  #version_label = 'V10_30Dec21'
  version_label = 'V10_30Dec21'
  dataset = 'D1'
  tagnano = '30Dec21'
  #tagnano = None
  tagflat = 'tag_fired_anyBParkHLT'
  #tagflat = 'fired_HLT_Mu9_IP6'

  TagAndProbeLauncher(is_data=is_data, out_label=out_label, version_label=version_label, ds=dataset, tagnano=tagnano, tagflat=tagflat).process()





