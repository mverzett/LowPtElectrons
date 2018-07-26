import pandas as pd
import numpy as np

cuts_df = pd.read_csv(
   '../../../RecoParticleFlow/PFTracking/data/Threshold.dat', 
   delim_whitespace=True, header=None
)

def baseline(df):
   for ibin in range(9):
      mask = (df.ibin == ibin)
      c_deta   = cuts_df.loc[ibin][0]
      c_dphi   = cuts_df.loc[ibin][1]
      c_ep     = cuts_df.loc[ibin][2]
      c_hits   = cuts_df.loc[ibin][3]
      c_chimin = cuts_df.loc[ibin][4]
      c_bdt    = cuts_df.loc[ibin][5]
      matching = (df[mask].trk_ecal_match > 0.5) &\
         (df[mask].trk_ecal_Deta < c_deta) & \
         (np.abs(df[mask].trk_ecal_Dphi) < c_dphi) & \
         (df[mask].e_over_p > c_ep) & \
         (df[mask].nhits > 10)
      togsf = (df[mask].nhits < c_hits) | \
         (df[mask].trk_chi2red > c_chimin)
      bdt_sel = np.invert(matching) & togsf & (df.bdtout > c_bdt)
      df.loc[mask, 'cutbased'] = (bdt_sel | matching)
      df.loc[mask, 'cutmatching'] = (matching)
      df.loc[mask, 'cutbdt'] = (bdt_sel)
      if 'bdt_thr' in df.columns:
         df.loc[mask, 'bdt_thr'] = c_bdt
