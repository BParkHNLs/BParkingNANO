import os
import sys
import ROOT

'''
  This script aims at retrieving what is the proportion of events at miniAOD level
  that are in the muon, electron channels respectively
'''

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to analyse gen particles from miniAOD samples', add_help=True)
  parser.add_argument('--f', type=str, dest='f', help='inputfile name', default='bparknano.root')
  return parser.parse_args()


def getChannelRate():

  opt = getOptions()

  #inputfile = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/signal_central/V02_generalStep/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL1p0_ctau1000p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/bparknano.root'
  inputfile = opt.f
  f = ROOT.TFile.Open(inputfile)
  tree = f.Get('Events')

  n_events = tree.GetEntries()
  print 'number of events ({}): {}'.format(inputfile, n_events)

  count_muon = 0
  count_electron = 0
  count_tot = 0

  for entry in tree:
    count_tot = count_tot + 1
    # searching for hnl
    has_hnl = False
    hnl_idx = -1
    b_idx = -1
    for icand in range(0, entry.nGenPart):
      if abs(entry.GenPart_pdgId[icand]) == 9900015:
        has_hnl = True
        hnl_idx = icand
        b_idx = entry.GenPart_genPartIdxMother[icand]
        break
    if not has_hnl: print 'WARNING: no hnl found'

    if has_hnl:
      # searching for hnl daughter
      daughter_ismuon = False
      daughter_iselectron = False
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == hnl_idx:
          if abs(entry.GenPart_pdgId[icand]) == 13:
            daughter_ismuon = True
          elif abs(entry.GenPart_pdgId[icand]) == 11: 
            daughter_iselectron = True

      # searching for hnl sister
      sister_ismuon = False
      sister_iselectron = False
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == b_idx:
          if abs(entry.GenPart_pdgId[icand]) == 13:
            sister_ismuon = True
          elif abs(entry.GenPart_pdgId[icand]) == 11: 
            sister_iselectron = True

      if daughter_ismuon and sister_ismuon: 
        count_muon = count_muon + 1
      elif (daughter_ismuon and sister_iselectron) or (daughter_iselectron and sister_ismuon):
        count_electron = count_electron + 1
      else:
        print 'Unknown channel'

  print 'filename: {}'.format(inputfile)
  print 'muon: {}'.format(round(float(count_muon) / float(count_tot), 2))
  print 'electron: {}'.format(round(float(count_electron) / float(count_tot), 2))
  print 'Check: muon + electron = {}'.format(float(count_muon) / float(count_tot) + float(count_electron) / float(count_tot))


def getEventsPerBspecies():

  #inputfile = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V44_Bu/mass2p0_ctau10p0/nanoFiles/merged/flat_bparknano_06Feb23.root'
  inputfile = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V42/mass4p0_ctau0p1/nanoFiles//merged/bparknano_generalstep.root'
  #inputfile = opt.f
  f = ROOT.TFile.Open(inputfile)
  tree = f.Get('Events')

  n_events = tree.GetEntries()
  print 'number of events ({}): {}'.format(inputfile, n_events)

  count_muon = 0
  count_electron = 0
  count_Bu_muon = 0
  count_Bd_muon = 0
  count_Bs_muon = 0
  count_Bu_electron = 0
  count_Bd_electron = 0
  count_Bs_electron = 0
  count_tot = 0

  for entry in tree:
    if count_tot%10000 == 0: print '--> processed {}% of the files'.format(round(float(count_tot)/n_events*100, 1))
    
    #TODO fetch number only per flavour channel
    count_tot += 1
    # searching for hnl
    has_hnl = False
    hnl_idx = -1
    b_idx = -1
    for icand in range(0, entry.nGenPart):
      if abs(entry.GenPart_pdgId[icand]) == 9900015:
        has_hnl = True
        hnl_idx = icand
        b_idx = entry.GenPart_genPartIdxMother[icand]
        break
    if not has_hnl: print 'WARNING: no hnl found'

    if has_hnl:
      # searching for hnl daughter
      daughter_ismuon = False
      daughter_iselectron = False
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == hnl_idx:
          if abs(entry.GenPart_pdgId[icand]) == 13:
            daughter_ismuon = True
          elif abs(entry.GenPart_pdgId[icand]) == 11: 
            daughter_iselectron = True

      # searching for hnl sister
      sister_ismuon = False
      sister_iselectron = False
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) in [11, 13] and entry.GenPart_genPartIdxMother[icand] == b_idx:
          if abs(entry.GenPart_pdgId[icand]) == 13:
            sister_ismuon = True
          elif abs(entry.GenPart_pdgId[icand]) == 11: 
            sister_iselectron = True

      if daughter_ismuon and sister_ismuon: 
        count_muon = count_muon + 1
        if abs(entry.GenPart_pdgId[b_idx]) == 521: count_Bu_muon += 1
        elif abs(entry.GenPart_pdgId[b_idx]) == 511: count_Bd_muon += 1
        elif abs(entry.GenPart_pdgId[b_idx]) == 531: count_Bs_muon += 1
      elif (daughter_ismuon and sister_iselectron) or (daughter_iselectron and sister_ismuon):
        count_electron = count_electron + 1
        if abs(entry.GenPart_pdgId[b_idx]) == 521: count_Bu_electron += 1
        elif abs(entry.GenPart_pdgId[b_idx]) == 511: count_Bd_electron += 1
        elif abs(entry.GenPart_pdgId[b_idx]) == 531: count_Bs_electron += 1
  
  print 'filename: {}'.format(inputfile)
  print 'muon rate: {}'.format(round(float(count_muon) / float(count_tot), 2))
  print 'electron rate: {}'.format(round(float(count_electron) / float(count_tot), 2))
  print 'number of Bu (muon): {}'.format(count_Bu_muon)
  print 'number of Bd (muon): {}'.format(count_Bd_muon)
  print 'number of Bs (muon): {}'.format(count_Bs_muon)
  print 'number of Bu (electron): {}'.format(count_Bu_electron)
  print 'number of Bd (electron): {}'.format(count_Bd_electron)
  print 'number of Bs (electron): {}'.format(count_Bs_electron)
  if count_Bu_muon + count_Bd_muon + count_Bs_muon + count_Bu_electron + count_Bd_electron + count_Bs_electron != count_tot: print '--> WARNING: numbers do not close'



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  getEventsPerBspecies()



