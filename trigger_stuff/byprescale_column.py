#! /bin/env python

import json
import argparse
import os
import sqlalchemy
import numpy as np
import cx_Oracle
import sys
import xml.etree.ElementTree as ET
import pandas as pd

class l1db :

    def __init__( self, user, password, database) :
        """
        database connection class
        """

        self.user = user
        self.password = password
        self.database = database
        self.connectionstring = '%s/%s@%s' % ( self.user, self.password, self.database )
        self.connection = cx_Oracle.connect(self.connectionstring)
        self.cursor = self.connection.cursor()

    def disconnect ( self ):
        """
        disconnection method
        """
        if self.connection is not None :
            self.connection.close()
        self.connection = None
        self.cursor = None

    def __del__(self):
        self.disconnect()

    def result ( self, columns ) :
        """
        fill map with query results
        """

        items = {}
        for item in columns :
            items[ item ] = []

        for row in self.cursor.fetchall() :
            for index, colname in enumerate( columns, start=0):
                items[ colname ].append( row[index] )

        return items


    def select( self, query ) :
        """
        execute query
        """
        try:
            print query
            self.cursor.execute(query)
        except cx_Oracle.DatabaseError, exc:
            error, = exc.args
            print >> sys.stderr, "Oracle-Error-Code:", error.code
            print >> sys.stderr, "Oracle-Error-Message:", error.message


    def fetchPrescales( self, runNumber ) :
        """
        select by timestamp query
        """
        columns = [ "LUMI_SECTION", "PRESCALE_INDEX", "PRESCALE_NAME" ]
        query = "select lumi.LUMI_SECTION, lumi.PRESCALE_INDEX, name.PRESCALE_NAME from CMS_UGT_MON.LUMI_SECTIONS lumi, CMS_UGT_MON.RUN_PRESCALENAMES name where lumi.PRESCALE_INDEX=name.PRESCALE_INDEX and lumi.RUN_NUMBER=name.RUN_NUMBER and lumi.RUN_NUMBER=%s order by lumi.LUMI_SECTION asc" % runNumber
        self.select(query)
        return self.result(columns)

    def fetch_fillinfo(self, run):
        columns = ['STARTTIME', 'STOPTIME', 'LHCFILL']
        query = 'select STARTTIME,STOPTIME,LHCFILL from CMS_WBM.RUNSUMMARY where RUNNUMBER=%s' % run
        self.select(query)
        return self.result(columns)

    def fetch_avgpu(self, run):
        columns = ['LSNUM', 'AVGPU']
        query = 'select LSNUM, AVGPU from cms_lumi_prod.online_result_6 where RUNNUM=%s' % run
        self.select(query)
        return self.result(columns)


class Ranges(object):
    def __init__(self, ranges):
        self.ranges = [(i[0], i[1]) for i in ranges]

    def __contains__(self, lumi):
        return any(i <= lumi <= j for i, j in self.ranges)

def call(cmd):
   print 'Executing: %s' % cmd
   retval = os.system(cmd)
   if retval != 0:
      raise RuntimeError('Command failed!')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('json', help='type of handle to be used')
args = parser.parse_args()

lumimap = json.load(open(args.json) )

#connect to trigger DB
engine = l1db( 'cms_trg_r', 'X3lmdvu4', 'cms_omds_lb' )
bril_engine = l1db('cms_runinfo_r', 'mickey2mouse', 'cms_orcon_adg')
data = pd.DataFrame(
    columns=[
        'fill',
        'fill_start',
        'fill_duration',
        'run',
        'prescale',
        'nlumis',
        'nlumiDCS',
        'date',
        'nPU', #o evtsize
        ]
    )

bycolumns = {}
def add_range(lumi, run, lrange):
   if lumi not in bycolumns: bycolumns[lumi] = {}
   if run not in bycolumns[lumi]: bycolumns[lumi][run] = []
   bycolumns[lumi][run].append(lrange)

remapping = {
   '1.0e34+BPH Parking' : '1.0e34',
   '1.4e34+Physics' : '1.4e34',
   '2.0e34+Physics' : '2.0e34',
   '1.4e34 backup'  : '1.4e34',
}

from pdb import set_trace
for run, lumis in lumimap.items():    
   irun = int(run)
   if irun < 315974: continue
   #get list of all prescales for such run invoking the dark forces (oracle)
   #suggested by our lord of DB and high sorcerer of the evil arts (Giuzz)
   result = engine.fetchPrescales(run)
   fillinfo = engine.fetch_fillinfo(run)
   puinfo = bril_engine.fetch_avgpu(run)
   puinfo = dict(zip(puinfo['LSNUM'], puinfo['AVGPU']))
   ## continue   
   lumirange = Ranges(lumis)
   full_lumis = result['PRESCALE_NAME']
   lumi_and_name = zip(result['LUMI_SECTION'], result['PRESCALE_NAME'])
   dcs_lumis = [j for i, j in lumi_and_name if i in lumirange]
   available_cols = set(full_lumis)
   for avail in available_cols:
       try:
           avg_pu = sum(
               puinfo[l] for l, n in lumi_and_name 
               if n == avail if l in puinfo
               )
           n_pu = len([
                   l for l, n in lumi_and_name 
                   if n == avail if l in puinfo
                   ])
       except:
           set_trace()
       fill_start = pd.to_datetime(fillinfo['STARTTIME'][0])
       fill_stop = pd.to_datetime(fillinfo['STOPTIME'][0])
       data.loc[data.shape[0]] = {
           'fill' : fillinfo['LHCFILL'][0],
           'fill_start' : fill_start,
           'fill_duration' : (fill_stop-fill_start).seconds,
           'run'  : irun,
           'prescale' : avail if avail not in remapping else remapping[avail],
           'nlumis'   : full_lumis.count(avail),
           'nlumiDCS' : dcs_lumis.count(avail),
           'date' : pd.to_datetime('2018/10/03'),
           'nPU'  : avg_pu/n_pu, #o evtsize
        }


bycolumn = pd.DataFrame(
    columns=[
        'prescale',
        'nlumis',
        'nlumiDCS',
        'nPU', #o evtsize
        ]
    )


for column, group in data.groupby('prescale'):
    bycolumn.loc[bycolumn.shape[0]] = {
        'prescale' : column,
        'nlumis'   : group.nlumis.sum(),
        'nlumiDCS' : group.nlumiDCS.sum(),
        'nPU'      : (group.nlumiDCS * group.nPU).sum()/group.nlumiDCS.sum(),
        }
bycolumn['evt_size'] = 200 + 20*bycolumn.nPU


proposal = { #rate, purity
   '2.2e34' : (0., 0.), 
   '2.0e34' : (1368, 0.63), 
   '1.6e34' : (1766, 0.72), 
   '1.4e34' : (2209, 0.76), 
   '1.2e34' : (2592, 0.79), 
   '1.0e34' : (2347, 0.78), 
   '8.0e33' : (3334*0.8/0.75, 0.57), 
   'HLTPhysics' : (0., 0.)#no parking 
}
lumisection_length = 23 #seconds
bycolumn['time'] = bycolumn.nlumiDCS*lumisection_length

hlt_options = pd.read_csv('machine_readable_selected.csv')
def get_rate_and_purity(col, seed, hlt):
    fseed = 'L1_SingleMu22 OR %s' % seed
    entry = hlt_options[
        (hlt_options.lumi == col) &
        (hlt_options.L1_seeds == fseed) &
        (hlt_options.HLT_Path == hlt)].iloc[0]
    return (entry.max_rate, entry.purity)

options = [
    ('current', {
            '2.2e34' : (0., 0.), 
            '2.0e34' : (0., 0.), 
            '1.6e34' : (0., 0.), 
            '1.4e34' : get_rate_and_purity(1.4, 'L1_SingleMu10er1p5', 'HLT_Mu9_IP6'),
            '1.2e34' : get_rate_and_purity(1.2, 'L1_SingleMu9er1p5' , 'HLT_Mu9_IP6'),
            '1.0e34' : get_rate_and_purity(1.0, 'L1_SingleMu8er1p5' , 'HLT_Mu9_IP6'),
            '8.0e33' : (3534*0.8/0.75, 0.57), #get_rate_and_purity(0.8, 'L1_SingleMu7er1p5' , 'HLT_Mu8_IP3'),
            'HLTPhysics' : (0., 0.)#no parking 
            }),
    ('+1.6', {
            '2.2e34' : (0., 0.), 
            '2.0e34' : (0., 0.), 
            '1.6e34' : get_rate_and_purity(1.6, 'L1_SingleMu12er1p5', 'HLT_Mu12_IP6'),
            '1.4e34' : get_rate_and_purity(1.4, 'L1_SingleMu10er1p5', 'HLT_Mu9_IP6'),
            '1.2e34' : get_rate_and_purity(1.2, 'L1_SingleMu9er1p5' , 'HLT_Mu9_IP6'),
            '1.0e34' : get_rate_and_purity(1.0, 'L1_SingleMu8er1p5' , 'HLT_Mu9_IP6'),
            '8.0e33' : (3534*0.8/0.75, 0.57),
            'HLTPhysics' : (0., 0.)#no parking 
            }),
    ('Mu12_IP6_all_the_way', {
            '2.2e34' : (0., 0.), 
            '2.0e34' : (0., 0.), 
            '1.6e34' : get_rate_and_purity(1.6, 'L1_SingleMu12er1p5', 'HLT_Mu12_IP6'),
            '1.4e34' : get_rate_and_purity(1.4, 'L1_SingleMu10er1p5', 'HLT_Mu12_IP6'),
            '1.2e34' : get_rate_and_purity(1.2, 'L1_SingleMu9er1p5' , 'HLT_Mu12_IP6'),
            '1.0e34' : get_rate_and_purity(1.0, 'L1_SingleMu8er1p5' , 'HLT_Mu12_IP6'),
            '8.0e33' : (3534*0.8/0.75, 0.57),
            'HLTPhysics' : (0., 0.)#no parking 
            }),
    ('+1.6_looser_1.0', {
            '2.2e34' : (0., 0.), 
            '2.0e34' : (0., 0.), 
            '1.6e34' : get_rate_and_purity(1.6, 'L1_SingleMu12er1p5', 'HLT_Mu12_IP6'),
            '1.4e34' : get_rate_and_purity(1.4, 'L1_SingleMu10er1p5', 'HLT_Mu9_IP6'),
            '1.2e34' : get_rate_and_purity(1.2, 'L1_SingleMu9er1p5' , 'HLT_Mu9_IP6'),
            '1.0e34' : get_rate_and_purity(1.0, 'L1_SingleMu8er1p5' , 'HLT_Mu9_IP5'),
            '8.0e33' : (3534*0.8/0.75, 0.57),
            'HLTPhysics' : (0., 0.)#no parking 
            }),
    ]

columns = [
 '2.2e34', 
 '2.0e34',
 '1.6e34', 
 '1.4e34', 
 '1.2e34', 
 '1.0e34',
 '8.0e33', 
 'HLTPhysics',
]
tot_lumis = bycolumn.nlumiDCS.sum()
M = lambda x: '%.1fM' % float(x*10**-6)
tot_exp_b = 0
tot_exp_t = 0
tot_opts = dict((i, {'b' : 0, 't' : 0}) for i, _ in options)
for col in columns:
   entry = bycolumn[bycolumn.prescale == col].iloc[0]
   nlumis = entry.nlumiDCS
   fraction = float(nlumis)/tot_lumis
   time = entry.time
   print 'Spent %.1f%% (%.2e s) in %s' % (100*fraction, time, col)
   rate, purity = proposal[col]
   exp_b = time * rate * purity
   exp_t = rate * entry.evt_size / 10**6
   tot_exp_b += exp_b
   tot_exp_t += float(exp_t*nlumis)
   print '   - #Bs expected: %s, throughput %.2f GB/s, rate %.2f kHz' % (
       M(exp_b), exp_t, rate/float(1000)
       )
   for name, rp_map in options:
       rate, purity = rp_map[col]
       b = time * rate * purity
       t = rate * entry.evt_size / 10**6
       tot_opts[name]['b'] += b
       tot_opts[name]['t'] += float(t*nlumis)
       print '   - %s: #Bs expected: %s (%+.2f%%), throughput %.2f GB/s (%+.2f%%), rate %.2f kHz' % (
           name, 
           M(b), (b/exp_b-1)*100 if exp_b else 0., 
           t, (t/exp_t-1)*100 if exp_t else 0.,
           rate/float(1000)
           )

   #print '#Bs expected: %s, recorded: %s (%+.2f%%)' % (M(exp), M(rec), (rec/exp-1)*100 if exp else 0)

print '-----------------------'
print 'Total:'
print '   - #Bs expected: %s, throughput (average) %.2f GB/s' % (
    M(tot_exp_b), tot_exp_t/tot_lumis
)
for name, _ in options:
    print '   - %s: #Bs expected: %s (%+.2f%%), throughput (average) %.2f GB/s (%+.2f%%)' % (
        name,
        M(tot_opts[name]['b']), (tot_opts[name]['b']/tot_exp_b-1)*100,
        tot_opts[name]['t']/tot_lumis, (tot_opts[name]['t']/tot_exp_t - 1)*100
    )

    
## effective_v1 = {
##    '2.2e34' : 0., 
##    '2.0e34' : 0.,
##    '1.6e34' : 0.,
##    '1.4e34' : 2592*1.4/1.2 * 0.79, 
##    '1.2e34' : 2592 * 0.79, 
##    '1.0e34' : 2300 * 0.78, 
##    '8.0e33' : 3534*0.8/0.75 * 0.57, #assume same purity as 1e34 for HLT_Mu9_IP6 
##    '1.0e34+BPH Parking' : 2300 * 0.78, 
##    '1.4e34+Physics' : 2592*1.4/1.2 * 0.79,
##    '2.0e34+Physics' : 0., 
##    '1.4e34 backup'  : 2592*1.4/1.2 * 0.79, 
##    'HLTPhysics' : 0., #no parking 
## }
## 
## nbs = 0
## nbs_effective = 0
## time_spent = {}
## for col, runs in bycolumns.iteritems():
##    pure_rate = proposal[col]
##    pure_effective_rate = effective_v1[col]
##    #print 'Using %f for %s' % (pure_rate, col)
##    ntime = 0
##    for lrange in [j for i in runs.values() for j in i]:
##       strt, stop = tuple(lrange)
##       length = ((stop+1) - strt)*lumisection_length
##       ntime += length
##       nbs += int(length*pure_rate) #convert to int to avoid problems
##       nbs_effective += int(length*pure_effective_rate)
##    time_spent[col] = ntime

