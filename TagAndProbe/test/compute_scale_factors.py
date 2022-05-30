import os
import sys
import glob
from os import path
import pandas as pd
from array import *
import ROOT



class ScaleFactorComputer(object):
  def __init__(self, data_label, mc_label, out_label):
    self.data_label = data_label
    self.mc_label = mc_label
    self.out_label = out_label


  def getDataSuffix(self, data_file):
    out_suffix = data_file[data_file.rfind(self.data_label)+len(self.data_label)+1:data_file.rfind('.root')]
    return out_suffix

  def getNEvents(self, out_suffix):
    return int(out_suffix[out_suffix.find('_n')+2:])


  def getItem(self, line, index):
    item = line[index:line.find(' ', index)]
    index = index + len(item)
    return item, index


  def computePerChunkSF(self, data_files, mc_file):
    '''
      Compute SF for each chunk of data separately
    '''

    print '\n --> will compute the SF per data chunk\n'

    scale_factors_perchunk = []
    n_events_perchunk = []
    for data_file in data_files:
      out_suffix = self.getDataSuffix(data_file) #data_file[data_file.rfind(self.data_label)+len(self.data_label)+1:data_file.rfind('.root')]

      if not path.exists('./results/{}/{}'.format(self.out_label, out_suffix)):
        os.system('mkdir -p ./results/{}/{}'.format(self.out_label, out_suffix))

      command_sf = 'sh submitter_sf.sh {data_file} {mc_file} {out_label}'.format(
          data_file = data_file,
          mc_file = mc_file, 
          out_label = self.out_label + '/' + out_suffix,
          )

      os.system(command_sf)

      scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_pt_eta_fit.txt'.format(self.out_label, out_suffix))
      #scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_pt_dxysig_fit.txt'.format(self.out_label, out_suffix))
      n_events_perchunk.append(self.getNEvents(out_suffix))

    return scale_factors_perchunk, n_events_perchunk


  def computeAverageScaleFactors(self, scale_factors_perchunk, n_events_perchunk, filename):
    '''
      The scale factors are defined as the weighted average of
      the scale factors in the different data chunks
    '''
    print '\n --> will compute the average SF\n'

    # create file containing the averaged scale factors
    average_sf_file = open(filename, 'w+')

    # get the total number of events
    n_events_tot = 0
    for ifile, sf_file in enumerate(scale_factors_perchunk):
      n_events_tot += n_events_perchunk[ifile]

    # compute the weighted scale factor per pt and eta bin
    content = []
    weighted_scale_factors = []
    for ifile, sf_file in enumerate(scale_factors_perchunk):

      # open per chunk sf file
      sf_file = open(sf_file)
      lines = sf_file.readlines()

      pt_bins = []
      eta_bins = []

      for line in lines:
        weighted_scale_factors_perchunk = {}
        ptmin, index = self.getItem(line, 0)
        ptmax, index = self.getItem(line, index+1)
        pt_bin = '{}_{}'.format(ptmin, ptmax)
        pt_bins.append(pt_bin)

        etamin, index = self.getItem(line, index+1)
        etamax, index = self.getItem(line, index+1)
        eta_bin = '{}_{}'.format(etamin, etamax)
        eta_bins.append(eta_bin)

        sf, index = self.getItem(line, index+1)
        weighted_sf = float(sf) * float(n_events_perchunk[ifile])
        weighted_scale_factors_perchunk['{}_{}'.format(pt_bin, eta_bin)] = weighted_sf
        weighted_scale_factors.append(weighted_scale_factors_perchunk)

        err, index = self.getItem(line, index+1)

    # compute the weighted average per pt and eta bin
    for pt_bin, eta_bin in zip(pt_bins, eta_bins):

      sum_weighted_scale_factor = 0.

      for weighted_scale_factor in weighted_scale_factors:
        if '{}_{}'.format(pt_bin, eta_bin) in weighted_scale_factor.keys():
          sum_weighted_scale_factor += weighted_scale_factor['{}_{}'.format(pt_bin, eta_bin)]

      average_scale_factor = sum_weighted_scale_factor / float(n_events_tot)
      pt_bin_min = pt_bin[:pt_bin.find('_')]
      pt_bin_max = pt_bin[pt_bin.find('_')+1:]
      eta_bin_min = eta_bin[:eta_bin.find('_')]
      eta_bin_max = eta_bin[eta_bin.find('_')+1:]
      content = '{} {} {} {} {}\n'.format(pt_bin_min, pt_bin_max, eta_bin_min, eta_bin_max, average_scale_factor) 
      #print content

      average_sf_file.write(content)

    average_sf_file.close()
    print '--> {} created'.format(filename)


  def process(self):
    # get data files
    data_files = [f for f in glob.glob('/work/anlyon/tag_and_probe/outfiles/{}/results*root'.format(self.data_label))]

    # get mc file
    mc_file = '/work/anlyon/tag_and_probe/outfiles/{lbl}/results_{lbl}_incl.root'.format(lbl=self.mc_label)

    # compute per chunk scale factors
    scale_factors_perchunk, n_events_perchunk = self.computePerChunkSF(data_files, mc_file)

    # compute the weighted average of the scale factors
    filename = './results/{}/scale_factors.txt'.format(self.out_label)
    self.computeAverageScaleFactors(scale_factors_perchunk, n_events_perchunk, filename)

    # create 2D histogram
    root_filename = './results/{}/scale_factors.root'.format(self.out_label)
    root_file = ROOT.TFile.Open(root_filename, "RECREATE")
    canv = ROOT.TCanvas('canv', 'canv', 800, 700)
    
    # open file with the scale factors
    sf_file = open(filename)
    lines = sf_file.readlines()

    pt_bins = [ 6. ] # set the lowest bin
    eta_bins = [ 0. ] # set the lowest bin

    scale_factor = {}

    for line in lines:
      ptmin, index = self.getItem(line, 0)
      ptmax, index = self.getItem(line, index+1)
      if float(ptmax) not in pt_bins: pt_bins.append(float(ptmax))

      etamin, index = self.getItem(line, index+1)
      etamax, index = self.getItem(line, index+1)
      if float(etamax) not in eta_bins: eta_bins.append(float(etamax))

      sf, index = self.getItem(line, index+1)
      scale_factor['{}_{}_{}_{}'.format(float(ptmin), float(ptmax), float(etamin), float(etamax))] = sf
      #print '{}_{}_{}_{}'.format(float(ptmin), float(ptmax), float(etamin), float(etamax))

    pt_bins = array('d', pt_bins)
    eta_bins = array('d', eta_bins)

    hist_scale_factor = ROOT.TH2D('hist_scale_factor', 'hist_scale_factor', len(pt_bins)-1, pt_bins, len(eta_bins)-1, eta_bins)

    for ipt, pt_bin in enumerate(pt_bins):
      for ieta, eta_bin in enumerate(eta_bins):
        if ipt == len(pt_bins)-1: continue
        if ieta == len(eta_bins)-1: continue
        #print '{} {} {}'.format(pt_bin, eta_bin, scale_factor['{}_{}_{}_{}'.format(pt_bin, pt_bins[ipt+1], eta_bin, eta_bins[ieta+1])])
        sf = float(scale_factor['{}_{}_{}_{}'.format(pt_bin, pt_bins[ipt+1], eta_bin, eta_bins[ieta+1])])
        hist_scale_factor.SetBinContent(ipt+1, ieta+1, sf)

    hist_scale_factor.SetOption("colztexte")
    hist_scale_factor.SetTitle("")
    hist_scale_factor.Write()
    hist_scale_factor.Draw()
    ROOT.gStyle.SetOptStat(0)
    
    canv.SaveAs('scale_factors.png')
    canv.SaveAs('scale_factors.pdf')
    root_file.Close()
    print '--> {} created'.format(root_filename)



if __name__ == "__main__":

  ROOT.gROOT.SetBatch(True)

  #data_label = 'scale_factors_tag_fired_anyBParkHLT_data_v2'
  #mc_label = 'scale_factors_tag_fired_anyBParkHLT_mc_v2'
  #out_label = 'scale_factors_tag_fired_anyBParkHLT_v2'

  #data_label = 'test_fullBPark_tag_fired_anyBParkHLT_max5e6'
  #mc_label = 'test_mc_tag_fired_anyBParkHLT'
  #out_label = 'test_fullBPark_tag_fired_anyBParkHLT_max5e6'

  data_label = 'test_D1_tag_fired_HLT_Mu12_IP6_pteta_max3e6'
  mc_label = 'test_mc_tag_fired_HLT_Mu12_IP6_pteta'
  out_label = 'test_D1_tag_fired_HLT_Mu12_IP6_pteta_max3e6'

  ScaleFactorComputer(data_label=data_label, mc_label=mc_label, out_label=out_label).process()


