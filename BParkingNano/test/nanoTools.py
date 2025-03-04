import os
import glob
import ROOT


class NanoTools(object):
  def getPointDirs(self, location):
    return [f for f in glob.glob(location+'/*')]


  def getNanoDirectories(self, location, prodlabel, dataset, dirtag=''):
    if dirtag == 'mcprivate':
      dirs = [f for f in glob.glob('{loc}/{pl}/*/nanoFiles'.format(loc=location, pl=prodlabel))]
    else:
      if dataset == None:
        dirs = [f for f in glob.glob('{loc}/{pl}/*'.format(loc=location, pl=prodlabel))]
      else:
        dirs = [f for f in glob.glob('{loc}/{pl}/{ds}'.format(loc=location, pl=prodlabel, ds=dataset))]
    if len(dirs) == 0:
      raise RuntimeError('No samples with the production label "{pl}" were found in {loc}'.format(pl=prodlabel if dataset==None else dataset+'_'+prodlabel, loc=location))
    return dirs


  def getLogDir(self, file_, prodlabel, tag, isData, dirlabel=''):
   #if isData: # probably to be modified for mc
   label = file_[file_.find('/',file_.find(prodlabel))+1:file_.find('Chunk')-1] 
   chunk = file_[file_.find('Chunk'):file_.find('bparknano')-1]
   if dirlabel == 'mcprivate':
     #logdir = './logs/{}/{}/{}'.format(prodlabel, label[0:label.find('/')], chunk) if tag == None else './logs/{}/{}_{}/{}'.format(prodlabel, label[0:label.find('/')], tag, chunk)
     #logdir = '/work/{}/logs/{}/{}/{}'.format(os.environ["USER"], prodlabel, label[0:label.find('/')], chunk) if tag == None else '/work/{}/logs/{}/{}_{}/{}'.format(os.environ["USER"], prodlabel, label[0:label.find('/')], tag, chunk)
     logdir = '/work/{}/CPV/logs/{}/{}/{}'.format(os.environ["USER"], prodlabel, label[0:label.find('/')], chunk) if tag == None else '/work/{}/CPV/logs/{}/{}_{}/{}'.format(os.environ["USER"], prodlabel, label[0:label.find('/')], tag, chunk)
   else:
     #logdir = './logs/{}/{}_{}/{}'.format(label, prodlabel, tag, chunk)
     logdir = '/work/{}/CPV/logs/{}/{}_{}/{}'.format(os.environ["USER"], label, prodlabel, tag, chunk) if tag != None else '/work/{}/CPV/logs/{}/{}/{}'.format(os.environ["USER"], label, prodlabel, chunk)
   return logdir


  def getFilesLocation(self, dirtag):  
    #location = '/pnfs/psi.ch/cms/trivcat/store/user/{usr}/BHNLsGen'.format(usr=os.environ["USER"])
    location = '/pnfs/psi.ch/cms/trivcat/store/user/{usr}/CPVGen'.format(usr=os.environ["USER"])
    if dirtag == 'data':
      location += '/data'
    elif dirtag == 'mccentral': 
      location += '/mc_central'
    elif dirtag == 'sigcentral': 
      location += '/signal_central'
    return location
  

  def getLogFile(self, logdir, file_):
    return '{}/nanostep_{}.log'.format(logdir, file_[file_.rfind('_nj')+1:file_.rfind('.root')])


  def getLocalMiniAODFiles(self, user, prodlabel, point, iscrab=False):
    if not iscrab:
      #pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/'.format(user, prodlabel, point)
      #pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/CPV/{}/{}/'.format(user, prodlabel, point)
      pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/CPVGen/{}/{}/'.format(user, prodlabel, point)
    else:
      #pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/*/*/*/'.format(user, prodlabel, point)
      #pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/CPV/{}/{}/*/*/*/'.format(user, prodlabel, point)
      pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/CPVGen/{}/{}/*/*/*/'.format(user, prodlabel, point)
    return [f for f in glob.glob(pointdir+'/step4_*.root')]


  def checkLocalFile(self, nanofile, cond=True, branch_check=False, branchname=''):
    rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
    if not rootFile: return False
    if cond and not rootFile.GetListOfKeys().Contains('Events'): return False
    else:
      if not branch_check: return True
      else:
        f = ROOT.TFile.Open(nanofile)
        tree = f.Get('Events')
        for branch in tree.GetListOfBranches():
          if(branch.GetName() == branchname):
            return True
        return False


  def checkFlatFile(self, flatfile, cond=True, branch_check=False, branchname=''):
    rootFile = ROOT.TNetXNGFile.Open(flatfile, 'r')
    if not rootFile: return False
    if cond and not rootFile.GetListOfKeys().Contains('signal_tree'): return False
    else:
      if not branch_check: return True
      else:
        f = ROOT.TFile.Open(flatfile)
        tree = f.Get('signal_tree')
        for branch in tree.GetListOfBranches():
          if(branch.GetName() == branchname):
            return True
        return False


  def checkFileExists(self, file_):
    import os.path
    from os import path
    return path.exists(file_)


  def getNFiles(self, file_):
    return sum(1 for line in open(file_))


  def getNExpectedNanoFiles(self, dir_):
    n_exp = dir_[dir_.rfind('_n')+2:len(dir_)]
    return int(n_exp)


  def getJobId(self, job):
    return int(job[job.find('job')+4:])


  def getJobIdsList(self, jobIds):
    listIds = ''
    for jobId in jobIds:
      listIds += '{}:'.format(jobId)
    return listIds[:len(listIds)-1]


  def getDataLabel(self, dataset):
    idx = dataset.find('-')
    return dataset.replace('/', '_')[1:idx]
  

  def getMCLabel(self, dataset):
    idx = dataset[1:].find('/')
    label = dataset[1:idx+1]
    if 'ext' in self.dataset: label += '_ext'
    return label


  def getStep(self, file_, iscrab=False): 
    if not iscrab:
      index = file_[file_.rfind('_nj')+3:file_.rfind('.root')]
    else:
      index = file_[file_.rfind('step4_')+6:file_.rfind('.root')]
    return index


  def getTag(self, tagnano, tagflat):
    if tagnano == None and tagflat == None:
      tag = 0
    elif tagnano != None and tagflat == None:
      tag = tagnano
    elif tagnano == None and tagflat != None:
      tag = tagflat
    else:
      tag = tagnano + '_' + tagflat
    return tag
    

