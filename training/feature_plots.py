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

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

class CMSJson(object):
    def __init__(self, fname):
        self.jmap = {int(i) : j for i,j in json.load(open('fill6371_JSON.txt')).iteritems()}
    
    def __contains__(self, key):
        run, lumi = key
        if run not in self.jmap: return False
        for lrange in self.jmap[run]:
            if lrange[0] <= lumi <= lrange[1]: return True
        return False
    
    def contains(self, run, lumi):
        cnt = lambda r, l: ((r, l) in self) #workaround to avoid self being picked-up by vectorize
        cnt = np.vectorize(cnt)
        return cnt(run, lumi)
    
    def __repr__(self):
        return 'CMSJson(%s)' % self.jmap

lumi_json = CMSJson('fill6371_JSON.txt')
data_signal = pd.DataFrame(root_numpy.root2array('data_features_signal.root'))
data_background = pd.DataFrame(root_numpy.root2array('data_features_background.root'))

print 'BEFORE JSON:'
print 'background:', data_background.shape
print 'signal:', data_signal.shape

lumi_mask = lumi_json.contains(data_signal.run, data_signal.lumi)
data_signal = data_signal[lumi_mask]
lumi_mask = lumi_json.contains(data_background.run, data_background.lumi)
data_background = data_background[lumi_mask]
data_signal['xy_sig'] = data_signal.dxy_err/data_signal.dxy
data_background['xy_sig'] = data_background.dxy_err/data_background.dxy

print 'AFTER JSON:'
print 'background:', data_background.shape
print 'signal:', data_signal.shape

plots = [u'nhits', u'ibin', u'high_purity',
       u'trk_ecal_match', u'bdtout', u'dxy', u'dxy_err', u'trk_pt', u'trk_inp',
       u'trk_outp', u'trk_eta', u'trk_ecal_Deta', u'trk_ecal_Dphi',
       u'e_over_p', u'trk_chi2red', u'gsf_success', u'gsf_dpt',
       u'trk_gsf_chiratio', u'gsf_chi2red', 'xy_sig']

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

for plot in plots:
    x_range = min(data_signal[plot].min(), data_background[plot].min()), \
        max(data_signal[plot].max(), data_background[plot].max())
    if plot in ranges: x_range = ranges[plot]
    plt.hist(data_signal[plot], bins=50, normed=True, histtype='step', label='electrons', range=x_range)
    plt.hist(data_background[plot], bins=50, normed=True, histtype='step', label='background', range=x_range)
    plt.legend(loc='best')
    plt.xlabel(plot if plot not in beauty else beauty[plot])
    plt.ylabel('A.U.')
    plt.show()
    plt.savefig('plots/%s.png' % plot)
    plt.savefig('plots/%s.pdf' % plot)
