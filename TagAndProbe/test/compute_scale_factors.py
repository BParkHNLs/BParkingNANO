import os
import sys
import glob
from os import path
import pandas as pd
from array import *
import ROOT



class ScaleFactorComputer(object):
  def __init__(self, data_label, mc_label, out_label, categorisation):
    self.data_label = data_label
    self.mc_label = mc_label
    self.out_label = out_label
    self.categorisation = categorisation
    self.eta_categories = ['0p00_0p50', '0p50_1p00', '1p00_1p50', '1p50_2p00']

    categories = ['pt_eta', 'pt_dxysig', 'pt_eta_dxysig']
    if self.categorisation not in categories:
      raise RuntimeError('Invalid categorisation. Please choose among {}'.format(categories))


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

      command_sf = 'sh submitter_sf.sh {data_file} {mc_file} {out_label} {cat}'.format(
          data_file = data_file,
          mc_file = mc_file, 
          out_label = self.out_label + '/' + out_suffix,
          cat = self.categorisation,
          )

      os.system(command_sf)

      if self.categorisation == 'pt_eta':
        scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_pt_eta_fit.txt'.format(self.out_label, out_suffix))
      elif self.categorisation == 'pt_dxysig':
        scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_pt_dxysig_fit.txt'.format(self.out_label, out_suffix))
      elif self.categorisation == 'pt_eta_dxysig':
        for eta_category in self.eta_categories:
          scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_eff_fit_eta_{}.txt'.format(self.out_label, out_suffix, eta_category))
          #scale_factors_perchunk.append('./results/{}/{}/scaleFactor_results_cat_pt_eta_dxysig_fit_eta_{}.txt'.format(self.out_label, out_suffix, eta_category))

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

    # compute the weighted scale factor per x and y bin
    content = []
    weighted_scale_factors = []
    for ifile, sf_file in enumerate(scale_factors_perchunk):

      # open per chunk sf file
      sf_file = open(sf_file)
      lines = sf_file.readlines()

      x_bins = []
      y_bins = []

      for line in lines:
        weighted_scale_factors_perchunk = {}
        xmin, index = self.getItem(line, 0)
        xmax, index = self.getItem(line, index+1)
        x_bin = '{}_{}'.format(xmin, xmax)
        x_bins.append(x_bin)

        ymin, index = self.getItem(line, index+1)
        ymax, index = self.getItem(line, index+1)
        y_bin = '{}_{}'.format(ymin, ymax)
        y_bins.append(y_bin)

        sf, index = self.getItem(line, index+1)
        weighted_sf = float(sf) * float(n_events_perchunk[ifile])
        weighted_scale_factors_perchunk['{}_{}'.format(x_bin, y_bin)] = weighted_sf
        weighted_scale_factors.append(weighted_scale_factors_perchunk)

        err, index = self.getItem(line, index+1)

    # compute the weighted average per x and y bin
    for x_bin, y_bin in zip(x_bins, y_bins):

      sum_weighted_scale_factor = 0.

      for weighted_scale_factor in weighted_scale_factors:
        if '{}_{}'.format(x_bin, y_bin) in weighted_scale_factor.keys():
          sum_weighted_scale_factor += weighted_scale_factor['{}_{}'.format(x_bin, y_bin)]

      average_scale_factor = sum_weighted_scale_factor / float(n_events_tot)
      x_bin_min = x_bin[:x_bin.find('_')]
      x_bin_max = x_bin[x_bin.find('_')+1:]
      y_bin_min = y_bin[:y_bin.find('_')]
      y_bin_max = y_bin[y_bin.find('_')+1:]
      content = '{} {} {} {} {}\n'.format(x_bin_min, x_bin_max, y_bin_min, y_bin_max, average_scale_factor) 
      #print content

      average_sf_file.write(content)

    average_sf_file.close()
    print '--> {} created'.format(filename)


  def createHistogram(self, sf_filename, name):
    # create 2D histogram
    root_file = ROOT.TFile.Open(name + '.root', "RECREATE")
    canv = ROOT.TCanvas('canv', 'canv', 800, 700)
    
    # open file with the scale factors
    sf_file = open(sf_filename)
    lines = sf_file.readlines()

    x_bins = [ 6. ] # set the lowest bin
    y_bins = [ 0. ] # set the lowest bin

    scale_factor = {}

    for line in lines:
      xmin, index = self.getItem(line, 0)
      xmax, index = self.getItem(line, index+1)
      if float(xmax) not in x_bins: x_bins.append(float(xmax))

      ymin, index = self.getItem(line, index+1)
      ymax, index = self.getItem(line, index+1)
      if float(ymax) not in y_bins: y_bins.append(float(ymax))

      sf, index = self.getItem(line, index+1)
      scale_factor['{}_{}_{}_{}'.format(float(xmin), float(xmax), float(ymin), float(ymax))] = sf
      #print '{}_{}_{}_{}'.format(float(xmin), float(xmax), float(ymin), float(ymax))

    x_bins = array('d', x_bins)
    y_bins = array('d', y_bins)

    hist_scale_factor = ROOT.TH2D('hist_scale_factor', 'hist_scale_factor', len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)

    for ix, x_bin in enumerate(x_bins):
      for iy, y_bin in enumerate(y_bins):
        if ix == len(x_bins)-1: continue
        if iy == len(y_bins)-1: continue
        #print '{} {} {}'.format(x_bin, y_bin, scale_factor['{}_{}_{}_{}'.format(x_bin, x_bins[ix+1], y_bin, y_bins[iy+1])])
        sf = float(scale_factor['{}_{}_{}_{}'.format(x_bin, x_bins[ix+1], y_bin, y_bins[iy+1])])
        hist_scale_factor.SetBinContent(ix+1, iy+1, sf)

    hist_scale_factor.SetOption("colztexte")
    hist_scale_factor.SetTitle("")
    hist_scale_factor.Write()
    hist_scale_factor.Draw()
    ROOT.gStyle.SetOptStat(0)
    
    canv.SaveAs(name + '.png')
    canv.SaveAs(name + '.pdf')
    root_file.Close()

    print '--> {}.root created'.format(name)


  def process(self):
    # get data files
    data_files = [f for f in glob.glob('/work/anlyon/tag_and_probe/outfiles/{}/results*root'.format(self.data_label))]

    # get mc file
    mc_file = '/work/anlyon/tag_and_probe/outfiles/{lbl}/results_{lbl}_incl.root'.format(lbl=self.mc_label)

    # compute per chunk scale factors
    scale_factors_perchunk, n_events_perchunk = self.computePerChunkSF(data_files, mc_file)

    if self.categorisation != 'pt_eta_dxysig':
      # compute the weighted average of the scale factors
      filename = './results/{}/scale_factors.txt'.format(self.out_label)
      self.computeAverageScaleFactors(scale_factors_perchunk, n_events_perchunk, filename)

      # create 2D histogram
      root_filename = './results/{}/scale_factors'.format(self.out_label)
      self.createHistogram(filename, root_filename)

    else:
      for eta_category in self.eta_categories:
        # compute the weighted average of the scale factors
        scale_factors = []
        for scale_factor in scale_factors_perchunk:
          if eta_category in scale_factor:
            scale_factors.append(scale_factor)
        filename = './results/{}/scale_factors_eta{}.txt'.format(self.out_label, eta_category)
        self.computeAverageScaleFactors(scale_factors, n_events_perchunk, filename)

        # create 2D histogram
        root_filename = './results/{}/scale_factors_eta{}'.format(self.out_label, eta_category)
        self.createHistogram(filename, root_filename)

    print '\nDone'


if __name__ == "__main__":

  ROOT.gROOT.SetBatch(True)

  data_label = 'test_fullBPark_tag_fired_anyBParkHLT_ptetadxysig_max5e6'
  mc_label = 'test_mc_tag_fired_anyBParkHLT_ptetadxysig'
  out_label = 'test_fullBPark_tag_fired_anyBParkHLT_ptetadxysig_max5e6'

  categorisation = 'pt_eta_dxysig'

  ScaleFactorComputer(data_label=data_label, mc_label=mc_label, out_label=out_label, categorisation=categorisation).process()


