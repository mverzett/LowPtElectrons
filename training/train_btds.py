import numpy as np
import matplotlib
matplotlib.use('Agg')
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
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
   'todo', nargs='+', choices=['reweight', 'train', 'evaluate', 'all'], 
   help='things to do (options: reweight, train, evaluate, all. If '
   'some task is not specified, the model is taken from the models directory'
)
parser.add_argument(
   '--mc', action='store_true', help='run training on MC'
)
args = parser.parse_args()

prefix = ('mc' if args.mc else 'data')
debug = False
data_signal = pd.DataFrame(
   root_numpy.root2array(
      'data/mc_features_signal.root' if args.mc else \
         'data/data_features_signal.root'
      )
)
data_background = pd.DataFrame(
   root_numpy.root2array(
      'data/mc_features_background.root' if args.mc else \
         'data/data_features_background.root'
      )
   )

bdt_features = [
   'nhits', 'high_purity', 'trk_pt', 'trk_inp',
   'trk_outp', 'trk_eta', 'trk_ecal_Deta', 'trk_ecal_Dphi',
   'e_over_p', 'trk_chi2red', 'gsf_dpt',
   'trk_gsf_chiratio', 'gsf_chi2red'
]

small_bdt_features = [
   'nhits', 'high_purity', 'trk_pt', 'trk_inp',
   'trk_outp', 'trk_eta', 'trk_ecal_Deta', 'trk_ecal_Dphi',
   'e_over_p', 'trk_chi2red', 
]

reweight_feats = ['trk_pt', 'trk_eta']

ranges = {
    'dxy' : (-0.1, 0.1),
    'dxy_err' : (0, 0.1),
    'trk_inp' : (0, 20),
    'trk_outp' : (0, 10),
    'trk_ecal_Deta' : (0, 0.2),
    'trk_ecal_Dphi' : (-0.2, 0.2),
    'e_over_p' : (0.1, 2.5),
    'trk_chi2red' : (0, 6),
    'gsf_dpt' : (0, 2),
    'trk_gsf_chiratio' : (0., 2),
    'gsf_chi2red' : (0, 4),
    'xy_sig' : (-6, 6),
    'nhits' : (0, 50),
}

beauty = {
    'dxy_err' : r'$\sigma$(dxy)', 
    'nhits' : r'\# of hits',
    'trk_pt' : r'p$_T$(ktf track)', 
    'log_trkpt' : r'log(p$_T$)(ktf track)', 
    'trk_inp' : r'p$_{in}$(ktf track)',
    'trk_outp' : r'p$_{out}$(ktf track)', 
    'trk_eta' : r'$\eta$(ktf track)', 
    'trk_ecal_Deta': '$\Delta\eta$(ECAL, ktf track)',
    'trk_ecal_Dphi' : '$\Delta\varphi$(ECAL, ktf track)',
    'e_over_p' : 'E/p', 
    'trk_chi2red' : '$\chi^2$(ktf track)/ndf', 
    'gsf_dpt' : r'p$_T$(gsf track)',
    'trk_gsf_chiratio' : '$\chi^2$(gsf track)/$\chi^2$(ktf track)', 
    'gsf_chi2red' : '$\chi^2$(gsf track)/ndf', 
    'xy_sig' : r'$\sigma$(dxy)/dxy',
}


#
# Merge samples
#
data_background['isE'] = 0
data_signal['isE'] = 1
data = pd.concat((data_background, data_signal))
data = data.sample(frac=1, random_state=42).reset_index(drop=True) #shuffle entries
data['log_trkpt'] = np.log(data.trk_pt)
#used later on but better having it here for data integrity
data['weight'] = 1. 
#sameprob = data_background.shape[0]/float(data_signal.shape[0])
#data.loc[(data.isE == 1), 'weight'] = sameprob
data['cutbased'] = 1. 
data['bigBDT'] = 1. 
data['smallBDT'] = 1. 

# Get train and test slice by hand to keep references and not copy datat
train = data[:int(.75*data.shape[0])]
test  = data[int(.75*data.shape[0]):]

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve, roc_auc_score
from hep_ml.reweight import GBReweighter
from sklearn.externals import joblib
import xgboost as xgb

#
# "Flatten" pt and eta
#
reweight_model_file = 'models/%s_reweighting.pkl' % prefix
subtrain = train[:int(.5*train.shape[0])]
subtest  = train[int(.5*train.shape[0]):]
if 'reweight' in args.todo or 'all' in args.todo:
   # Start checking the separation at start
   pre_separation = GradientBoostingClassifier(
      n_estimators=1 if debug else 50, 
      max_depth=4, random_state=42, verbose=1
   )
   #set_trace()
   pre_separation.fit(subtrain[reweight_feats], subtrain[['isE']], sample_weight=subtrain.weight)
   test_proba = pre_separation.predict_proba(subtest[reweight_feats])[:, 1]
   roc_pre = roc_curve(subtest[['isE']],  test_proba, sample_weight=subtest.weight)[:2]
   auc_pre = roc_auc_score(subtest[['isE']],  test_proba, sample_weight=subtest.weight)
   
   #run reweighting -- not working on MC for some reason
   reweighter = GBReweighter(
      n_estimators=1 if debug else 30, 
      max_depth=4, learning_rate=0.1
      )
   reweighter.fit(subtrain[subtrain.isE == 1][reweight_feats], subtrain[subtrain.isE == 0][reweight_feats]) #make electrons look like tracks
   
   #run weights FOR EVERYTHING!
   for df in [data, subtrain, subtest]: 
      weights = reweighter.predict_weights(df[df.isE == 1][reweight_feats]) #1/w to be used
      df.loc[df.isE == 1, 'weight'] = weights

   #save reweighter
   joblib.dump(reweighter, reweight_model_file, compress=True)
   
   # Check that sepratation vanishes
   post_separation = GradientBoostingClassifier(
      n_estimators=1 if debug else 50, 
      max_depth=4, random_state=42, verbose=1)
   post_separation.fit(
      subtrain[reweight_feats], 
      subtrain[['isE']], 
      sample_weight=subtrain.weight
   )
   test_proba = post_separation.predict_proba(subtest[reweight_feats])[:, 1]
   roc_post = roc_curve(subtest[['isE']],  test_proba, sample_weight=subtest.weight)[:2]
   auc_post = roc_auc_score(subtest[['isE']],  test_proba, sample_weight=subtest.weight)
   
   # make plots
   plt.figure(figsize=[8, 8])
   plt.plot(*roc_pre, label='Initial separation (%.3f)' % auc_pre)
   plt.plot(*roc_post, label='Separation after reweighting (%.3f)' % auc_post)
   plt.plot([0, 1], [0, 1], 'k--')
   plt.xlabel('FPR')
   plt.ylabel('TPR')
   plt.legend(loc='best')
   plt.plot()
   plt.savefig('plots/%s_reweighting.png' % prefix)
   plt.savefig('plots/%s_reweighting.pdf' % prefix)
   plt.clf()
   
   for plot in reweight_feats:
      x_range = min(data_signal[plot].min(), data_background[plot].min()), \
         max(data_signal[plot].max(), data_background[plot].max())
      if plot in ranges: x_range = ranges[plot]
      plt.hist(
         data[data.isE == 1][plot], bins=50, normed=True, 
         histtype='step', label='electrons', range=x_range, weights=data[data.isE == 1].weight
         )
      plt.hist(
         data[data.isE == 0][plot], bins=50, normed=True, 
         histtype='step', label='background', range=x_range, weights=data[data.isE == 0].weight
         )
      plt.legend(loc='best')
      plt.xlabel(plot if plot not in beauty else beauty[plot])
      plt.ylabel('A.U.')   
      plt.savefig('plots/%s_reweight_%s.png' % (prefix, plot))
      plt.savefig('plots/%s_reweight_%s.pdf' % (prefix, plot))
      plt.clf()
   if (auc_pre - auc_post) < 0.05:
      print 'Something went wrong with the reweighting, terminating here! (%.3f --> %.3f)' % (auc_pre, auc_post)
      exit()
else:
   reweighter = joblib.load(reweight_model_file)

#run weights FOR EVERYTHING!
for df in [data, train, test]: #ugly, but waiting for a fix
   weights = reweighter.predict_weights(df[df.isE == 1][reweight_feats]) #1/w to be used
   df.loc[df.isE == 1, 'weight'] = weights
print 'Reweighting done! (%.3f --> %.3f)' % (auc_pre, auc_post)

#
# Train BDTs
#
full_model = 'models/%s_full_BDT.pkl' % prefix
partial_model = 'models/%s_partial_BDT.pkl' % prefix
if 'train' in args.todo or 'all' in args.todo:   
   big_clf = xgb.XGBClassifier(
      #basic stuff
      max_depth=3, learning_rate=0.1, n_estimators=1 if debug else 100, 
      #many different ways of regularization
      gamma=0, min_child_weight=1, max_delta_step=0, subsample=1, 
      colsample_bytree=1, colsample_bylevel=1, 
      reg_alpha=0, reg_lambda=1, 
      #running settings and weight balancing
      silent=False, nthread=6, scale_pos_weight=1, 
   )
   big_clf.fit(train[bdt_features], train.isE, sample_weight=train.weight)
   joblib.dump(big_clf, full_model, compress=True)
   
   small_clf = xgb.XGBClassifier(
      #basic stuff
      max_depth=3, learning_rate=0.1, n_estimators=1 if debug else 100, 
      #many different ways of regularization
      gamma=0, min_child_weight=1, max_delta_step=0, subsample=1, 
      colsample_bytree=1, colsample_bylevel=1, 
      reg_alpha=0, reg_lambda=1, 
      #running settings and weight balancing
      silent=False, nthread=6, scale_pos_weight=1, 
   )
   small_clf.fit(train[small_bdt_features], train.isE, sample_weight=train.weight)
   joblib.dump(small_clf, partial_model, compress=True)
else:
   big_clf = joblib.load(full_model)
   small_clf = joblib.load(partial_model)

test['bigBDT'] = big_clf.predict_proba(test[bdt_features])[:, 1]
test['smallBDT'] = small_clf.predict_proba(test[small_bdt_features])[:, 1]
print 'Training done!'

#
# plot performance
#
if 'evaluate' in args.todo or 'all' in args.todo:
   compare_signal = pd.DataFrame(
      root_numpy.root2array(
         'data/mc_features_signal.root' if not args.mc else \
            'data/data_features_signal.root'
         )
   )
   compare_background = pd.DataFrame(
      root_numpy.root2array(
         'data/mc_features_background.root' if not args.mc else \
            'data/data_features_background.root'
         )
      )
   compare_background['isE'] = 0
   compare_signal['isE'] = 1
   compare = pd.concat((compare_background, compare_signal))
   compare = compare.sample(frac=1, random_state=42).reset_index(drop=True) #shuffle entries
   compare['bigBDT']   = big_clf.predict_proba(compare[bdt_features])[:, 1] 
   compare['smallBDT'] = small_clf.predict_proba(compare[small_bdt_features])[:, 1] 
   
   ## cuts = pd.read_csv(
   ##    '../../../RecoParticleFlow/PFTracking/data/Threshold.dat', 
   ##    delim_whitespace=True, header=None
   ## )
   ## for ibin in range(9):
   ##    mask = (test.ibin == ibin)
   ##    c_deta = cuts[ibin][0]
   ##    c_dphi = cuts[ibin][1]
   ##    c_ep = cuts[ibin][2]
   ##    test.loc[mask, 'cutbased'] = (test[mask].trk_ecal_Deta < c_deta) & \
   ##       (np.abs(test[mask].trk_ecal_Dphi) < c_dphi) & \
   ##       (test[mask].e_over_p > c_ep) & \
   ##       (test[mask].nhits > 10)
   ## cut_based = (
   ##    [float((test.cutbased > 0.5).sum())/test[test.isE < 0.5].shape[0]],
   ##    [float((test.cutbased > 0.5).sum())/test[test.isE > 0.5].shape[0]],
   ##    )
   
   full_roc = roc_curve(test.isE, test.bigBDT)[:2]
   small_roc = roc_curve(test.isE, test.smallBDT)[:2]
   old_roc = roc_curve(test.isE, test.bdtout)[:2]
   
   compare_full_roc  = roc_curve(compare.isE, compare.bigBDT)[:2]
   compare_small_roc = roc_curve(compare.isE, compare.smallBDT)[:2]
   compare_old_roc   = roc_curve(compare.isE, compare.bdtout)[:2]
   other_label = 'data' if args.mc else 'mc'
   
   # make plots
   plt.figure(figsize=[8, 8])
   plt.title('%s training' % prefix)
   plt.plot(*full_roc , color='r', label='Retraining, full BDT (%s)'   % prefix)
   plt.plot(*old_roc  , color='g', label='Old training, full BDT (%s)' % prefix)
   plt.plot(*small_roc, color='b', label='Retraining, small BDT (%s)'  % prefix)
   plt.plot(*compare_full_roc , color='r', linestyle='--', label='Retraining, full BDT (%s)'   % other_label)
   plt.plot(*compare_old_roc  , color='g', linestyle='--', label='Old training, full BDT (%s)' % other_label)
   plt.plot(*compare_small_roc, color='b', linestyle='--', label='Retraining, small BDT (%s)'  % other_label)
   plt.plot([0, 1], [0, 1], 'k--')
   plt.xlabel('FPR')
   plt.ylabel('TPR')
   plt.legend(loc='best')
   plt.plot()
   plt.savefig('plots/%s_training.png' % prefix)
   plt.savefig('plots/%s_training.pdf' % prefix)
   plt.clf()
