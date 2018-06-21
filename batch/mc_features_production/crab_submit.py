from WMCore.Configuration import Configuration

config = Configuration()
config.section_('General')
config.General.requestName = '2018Jun19_LowPtElectrons'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/m/mverzett/RK94New/src/CMGTools/LowPtElectrons/run/raw2features.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['outname=mc_features.root', 'data=RAWMCTest', 'globalTag=94X_mc2017_realistic_v12', 'isMC=True']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['mc_features_background.root', 'mc_features_signal.root']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/cmst3/user/mverzett/'
config.Data.splitting = 'FileBased'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.Data.inputDataset = '/BToKee_Pythia/tstreble-BToKee_Pythia_PUMix_18_03_17-c9b9e020b5bce5ee6bee9ef5f38c415a/USER'
config.Site.whitelist = ['T2_UK_London_IC']
config.Data.inputDBS = 'phys03'
