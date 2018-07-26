import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
from cmsjson import CMSJson

import matplotlib.pyplot as plt
import ROOT
import root_numpy
import rootpy
import json
import pandas as pd
from matplotlib import rc
from pdb import set_trace
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from baseline import baseline

debug = False
print 'Getting the data'
input_files = {
   'mc' : {
      'sig' : 'data/mc_features_signal_v3.root',
      'bkg' : 'data/mc_features_background_v3.root',
      },
   'data' : {
      'sig' : 'data/data_features_signal_v2.root',
      'bkg' : 'data/data_features_background_v2.root',
      },
}

data_signal = pd.DataFrame(
   root_numpy.root2array(
      input_files['data']['sig'],
      stop = 1000
      )
   )
data_background = pd.DataFrame(
   root_numpy.root2array(
      input_files['data']['bkg'],
      stop = 1000
      )
   )
mc_signal = pd.DataFrame(
   root_numpy.root2array(
      input_files['mc']['sig'],
      stop = 1000
      )
   )
mc_background = pd.DataFrame(
   root_numpy.root2array(
      input_files['mc']['bkg'],
      stop = 1000
      )
   )

missing = set(list(mc_signal.columns)) - set(list(mc_background.columns)) 
for i in missing:
   mc_background[i] = 0

lumi_json = CMSJson('fill6371_JSON.txt')
lumi_mask = lumi_json.contains(data_signal.run, data_signal.lumi)
data_signal = data_signal[lumi_mask]
lumi_mask = lumi_json.contains(data_background.run, data_background.lumi)
data_background = data_background[lumi_mask]
data_background['isE'] = 0
data_signal['isE'] = 1
data = pd.concat((data_background, data_signal))
data = data.sample(frac=1, random_state=42).reset_index(drop=True) #shuffle entries
#used later on but better having it here for data integrity
#sameprob = data_background.shape[0]/float(data_signal.shape[0])
#data.loc[(data.isE == 1), 'weight'] = sameprob
data['cutbased'] = False
data['cutmatching'] = False
data['cutbdt'] = False
baseline(data)


mc_background['isE'] = 0
mc_signal['isE'] = 1
mc = pd.concat((mc_background, mc_signal))
mc = mc.sample(frac=1, random_state=42).reset_index(drop=True) #shuffle entries
mc['cutbased'] = False
mc['cutmatching'] = False
mc['cutbdt'] = False
baseline(mc)
X_ = lambda x: [i for i, _, _ in x]
Y_ = lambda x: [i for _, i, _ in x]
Z_ = lambda x: [i for _, _, i in x]

for var, binning, xlegend in [
   ('trk_pt', np.arange(1, 11, 1), 'ktf track pT'),
   ('trk_eta', np.arange(-2.5, 2.6, 0.5), 'ktf track eta')]:
   lows = binning[:-1]
   highs = binning[1:]
   cutb  = []
   comp_cutb  = []
   for low, high in zip(lows, highs):
      mask = (low <= data[var]) & (data[var] < high)
      masked = data[mask]
      eff = float(masked[masked.isE > 0.5].cutbased.sum())/masked[masked.isE > 0.5].shape[0]
      ntot = masked[masked.isE > 0.5].shape[0]
      cutb.append((
            (low+high)/2.,
            eff,
            np.sqrt(eff * (1- eff) / ntot)
            ))
      
      mask = (low <= mc[var]) & (mc[var] < high)
      masked = mc[mask]
      eff = float(masked[masked.isE > 0.5].cutbased.sum())/masked[masked.isE > 0.5].shape[0]
      ntot = masked[masked.isE > 0.5].shape[0]
      comp_cutb.append((
            (low+high)/2.,
            eff,
            np.sqrt(eff * (1- eff) / ntot)
            ))
   
   plt.clf()
   plt.errorbar(X_(cutb) , Y_(cutb) , yerr=Z_(cutb) , color='g', fmt='-o', label='Current baseline (data)')
   plt.errorbar(X_(comp_cutb) , Y_(comp_cutb) , yerr=Z_(comp_cutb) , color='g', fmt='--o', label='Current baseline (mc)')
   plt.xlabel(xlegend)
   plt.ylabel('Efficiency')
   plt.legend(loc='best')
   plt.grid(which='both')
   plt.savefig('plots/baseline_vs_%s.png' % (var))
   plt.savefig('plots/baseline_vs_%s.pdf' % (var))
   plt.clf()



#TODO
gen_electrons = mc[mc.isE > 0.5]
for var, binning, xlegend in [
   ('gen_pt', np.arange(1, 11, 1), 'gen electron pT')]:
   lows = binning[:-1]
   highs = binning[1:]
   cutb  = []
   comp_cutb  = []
   for low, high in zip(lows, highs):
      mask = (low <= gen_electrons[var]) & (gen_electrons[var] < high)
      masked = gen_electrons[mask]
      eff = float(  masked.cutbased.sum())/masked.shape[0]
      ntot = masked.shape[0]
      comp_cutb.append((
            (low+high)/2.,
            eff,
            np.sqrt(eff * (1- eff) / ntot)
            ))
   
   plt.clf()
   plt.errorbar(X_(comp_cutb) , Y_(comp_cutb) , yerr=Z_(comp_cutb) , color='g', fmt='--o', label='Current baseline (mc)')
   plt.xlabel(xlegend)
   plt.ylabel('Efficiency')
   plt.grid(which='both')
   plt.legend(loc='best')
   plt.savefig('plots/baseline_vs_%s.png' % (var))
   plt.savefig('plots/baseline_vs_%s.pdf' % (var))
   plt.clf()

plt.scatter(
   gen_electrons.gen_pt, gen_electrons.trk_pt, 
   s=np.ones(gen_electrons.gen_pt.shape[0])
)
plt.grid(which='both')
plt.xlabel('gen pT')
plt.ylabel('ktf track pT')
plt.savefig('plots/trk_pT_vs_gen_pT.png')
plt.savefig('plots/trk_pT_vs_gen_pT.pdf')
plt.clf()

   
