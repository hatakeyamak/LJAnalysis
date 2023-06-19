#! /usr/bin/env python

## Plots for lepton-quark scattering to lepton + jet + no MET final state

import os
import sys
import pwd
import math
import time
import subprocess
import numpy as np
import ROOT as R

## For filling trees
from array import array
from collections import OrderedDict

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## Load local python modules
sys.path.insert(0, '%s/python' % os.getcwd())
import EventWeights
EVT_WGT = EventWeights.GetWeights()
import EventSelection
EVT_SEL = EventSelection.GetSelection()

## User configuration
TREE_OUT = False

PRT_EVT = 1  ## Print every Nth event
MAX_EVT = 1000      ## Number of events to process
# LUMI    = 59830  ## 2018 integrated lumi (pb-1), certified "Good"
LUMI    = 6900  ## 2018C integrated lumi (pb-1), certified "Good"

LABEL = 'XXX'
# CAT   = 'SingleMu'
# CAT   = 'SingleMu_oneLep'
# CAT   = 'SingleMu_oneLep_lepMVA0p8'
# CAT   = 'SingleMu_oneLep_lepMVA0p8_pMET30'
CAT   = 'SingleMu_oneLep_lepMVA0p8_pMET30_ge1j'
# CAT   = 'SingleMu_oneLep_lepMVA0p8_pMET30_ge1j_cenJet30'
NJOB  = 1
IJOB  = 0
USER  = pwd.getpwuid(os.getuid())[0]


if len(sys.argv) > 1:
    print '\nLABEL changed from %s to %s' % (LABEL, str(sys.argv[1]))
    LABEL = str(sys.argv[1])
if len(sys.argv) > 2:
    print '\nCAT changed from %s to %s' % (CAT, str(sys.argv[2]))
    CAT = str(sys.argv[2])
if len(sys.argv) > 3:
    print '\nMAX_EVT changed from %d to %d' % (MAX_EVT, int(sys.argv[3]))
    MAX_EVT = int(sys.argv[3])
    PRT_EVT = MAX_EVT / 100
if len(sys.argv) > 4:
    print '\nWill split into %d jobs, run job #%d' % (int(sys.argv[4]), int(sys.argv[5]))
    NJOB = int(sys.argv[4])
    IJOB = int(sys.argv[5])

if TREE_OUT:
    print '\n\n*** WARNING!!! Only TTrees will be output, not histograms.  Quit now if you want histograms! ***\n\n'
    time.sleep(5)

    
def main():

    print '\nInside LeptonJet\n'

    in_file_names = []
    in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/'
    if LABEL == 'SingleMuon_2018C':       in_dir += 'data/SingleMuon/UL2018_MiniAODv2_NanoAODv9-v2/2018C/'
    if LABEL == 'QCD_MuEnriched':         in_dir += 'MC/QCD/QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8/UL18NanoAODv9/'
    if LABEL == 'WJetsToLNu':             in_dir += 'MC/WJets/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/UL18NanoAODv9/'
    if LABEL == 'DYJetsToLL_M-10to50':    in_dir += 'MC/ZJets/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/UL18NanoAODv9/'
    if LABEL == 'DYJetsToLL_M-50':        in_dir += 'MC/ZJets/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/UL18NanoAODv9/'
    if LABEL == 'LQsc-SM_ktMin_20':       in_dir += 'MC/LeptonJet/LQsc-SM_ktMin_20/UL18NanoAODv9/split/'
    if LABEL == 'LQsc-SM_up_mu_ktMin_20': in_dir += 'MC/LeptonJet/LQsc-SM_up_mu_ktMin_20/UL18NanoAODv9/'
    if LABEL == 'LQsc-SM_up_mu_ktMin_60': in_dir += 'MC/LeptonJet/LQsc-SM_up_mu_ktMin_60/UL18NanoAODv9/'

    if not LABEL.split('_')[0] in in_dir:
        print '\n\n***  TRIED TO APPLY LABEL %s TO DIRECTORY %s  ***' % (LABEL, in_dir)
        print '***  CORRECT? QUITTING!  ***'
        sys.exit()

    nFiles = len(subprocess.check_output(['ls', in_dir]).splitlines()) - 1
    nJobs  = min(NJOB, nFiles)
    print '\n\n*** Found %d total input files ***\n\n' % nFiles
    if IJOB >= nJobs:
        print '\n\n***  TRIED TO RUN IJOB = %d WITH ONLY %d JOBS TOTAL!!! ***' % (IJOB, nJobs)
        sys.exit()
    iFile = 0
    for f_name in subprocess.check_output(['ls', in_dir]).splitlines():
        if not '.root' in f_name: continue
        # if not 'D7033346-A91B-5440-9B0A-8FCB0D90EF4B.root' in f_name: continue  ## 3M evt
        # if not 'D0D7BD3D-5626-6C4B-B9F4-8DBD9AE223A9.root' in f_name: continue  ## 350k evt
        # if not '2B07B4C0-852B-9B4F-83FA-CA6B047542D1.root' in f_name: continue  ## 20k evt
        if 'LQsc-SM_up_mu' in f_name and not 'GTfix' in f_name: continue
        if (iFile % nJobs) == IJOB:
            in_file_names.append(in_dir+f_name)
            print 'Appending file: %s' % in_file_names[-1]
        iFile += 1

    out_dir = 'plots/'
    if TREE_OUT: out_dir = '/cms/data/store/user/%s/Trees/LeptonJet/' % USER
    if nJobs > 1:
        out_dir += 'split/'
    out_file_str = out_dir+'LeptonJet_%s_%s' % (CAT, LABEL)
    if MAX_EVT > 0:
        if MAX_EVT >= 1000000: out_file_str += '_%dM' % (MAX_EVT / 1000000)
        else:                  out_file_str += '_%dk' % (MAX_EVT / 1000)
    if nJobs > 1:
        out_file_str += '_%d_%03d' % (nJobs, IJOB)
    out_file_str += '_TEST'
    out_file_str += '.root'
    out_file = R.TFile(out_file_str,'recreate')

    if TREE_OUT:
        out_tree = R.TTree("vars", "vars")
        int_vars = OrderedDict()
        flt_vars = OrderedDict()
        ints = ['proc']
        flts = ['wgt',
                'CoM_log2_pt_lj','CoM_log2_pz_lj','CoM_log2_p_jet','CoM_dTheta_pz_jet','CoM_dPhi_pt_jet',
                'Det_log2_pt_lep','Det_log2_pt_jet','Det_log2_pt_pMET','Det_log2_ptRel_lj','Det_log2_ptRel_lv',
                'dPhi_lv','dPhi_lj','dPhi_jv','dPhi_ljv']

        for var in ints:
            int_vars[var] = array('i', [0])
            out_tree.Branch(var, int_vars[var], var+'/I')
        for var in flts:
            flt_vars[var] = array('f', [0.])
            out_tree.Branch(var, flt_vars[var], var+'/F')


    isData = False
    if 'SingleMuon' in in_dir or 'EGamma' in in_dir: isData = True

    chains = {}
    chains['Events'] = 0

    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )
            print 'Added TChain %s' % key


    ## Histogram binning for each variable
    bins = {}
    bins['pt']     = [200,  0, 1000]
    bins['pz']    = [2000, -5000, 5000]
    bins['eta']    = [200, -5,    5]
    bins['phi']    = [128, -3.2, 3.2]
    bins['charge'] = [  5, -2.5, 2.5]
    bins['flavor'] = [  5, -2.5, 2.5]
    bins['sum']    = [100,  0, 1000]
    bins['dR']     = [100,  0,   10]
    bins['dEta']   = [100,  0,   10]
    bins['dPhi']   = [ 32,  0,  3.2]
    bins['mass']   = [200,  0, 2000]
    bins['MT']     = [200,  0, 2000]
    bins['HT']     = [200,  0, 2000]


    ## Sets of objects to be plotted
    Objs = ['lep1','lep2','MET','pMET','lv','lj','jv','ll','jj','ljj','ljv','jet1','jet2','jet3','evt']
    Vars = ['pt','pz','eta','phi','charge','flavor','dR','dEta','dPhi','mass','MT','HT']

    ## Book histograms (most will not be used, and will be deleted at the end)
    hst = {}
    ## All combinations of one and two objects for each variable
    for var in Vars:
        for obj1 in Objs:
            hst['%s_%s' % (var, obj1)] = R.TH1D( 'h_%s_%s' % (var, obj1), '%s %s' % (obj1, var),
                                                bins[var][0], bins[var][1], bins[var][2] )
            for obj2 in Objs:
                hst['%s_%s_%s' % (var, obj1, obj2)] = R.TH1D( 'h_%s_%s_%s' % (var, obj1, obj2),
                                                              '%s(%s, %s)' % (var, obj1, obj2),
                                                              bins[var][0], bins[var][1], bins[var][2] )


    hst['nLep'] = R.TH1D('h_nLep', 'h_nLep',  5, -0.5, 4.5)
    hst['nMu']  = R.TH1D('h_nMu',  'h_nMu',   5, -0.5, 4.5)
    hst['nEle'] = R.TH1D('h_nEle', 'h_nEle',  5, -0.5, 4.5)
    hst['nJet'] = R.TH1D('h_nJet', 'h_nJet', 10, -0.5, 9.5)

    del hst['mass_ll']
    hst['mass_ll'] = R.TH1D('h_mass_ll', 'h_mass_ll', 200, 0, 400)

    hst['charge_sum_ll'] = R.TH1D('h_charge_sum_ll', 'h_charge_sum_ll', 5, -2.5, 2.5)
    hst['flavor_sum_ll'] = R.TH1D('h_flavor_sum_ll', 'h_flavor_sum_ll', 5, -2.5, 2.5)

    hst['mvaTTH_lep1']   = R.TH1D('h_mvaTTH_lep1',   'h_mvaTTH_lep1',   40, -1, 1)
    hst['mvaTTH_lep2']   = R.TH1D('h_mvaTTH_lep2',   'h_mvaTTH_lep2',   40, -1, 1)
    hst['mvaTTH_min_ll'] = R.TH1D('h_mvaTTH_min_ll', 'h_mvaTTH_min_ll', 40, -1, 1)

    hst['miniIso_lep1']   = R.TH1D('h_miniIso_lep1',   'h_miniIso_lep1',   11, 0, 0.22)
    hst['miniIso_lep2']   = R.TH1D('h_miniIso_lep2',   'h_miniIso_lep1',   11, 0, 0.22)
    hst['miniIso_max_ll'] = R.TH1D('h_miniIso_max_ll', 'h_miniIso_max_ll', 16, 0, 0.8)

    hst['ID_lep1'] = R.TH1D('h_ID_lep1', 'h_ID_lep1', 4, -0.5, 3.5)
    hst['ID_lep2'] = R.TH1D('h_ID_lep2', 'h_ID_lep2', 4, -0.5, 3.5)

    hst['SIP_lep1']   = R.TH1D('h_SIP_lep1',   'h_SIP_lep1',   16, 0, 8.0)
    hst['SIP_lep2']   = R.TH1D('h_SIP_lep2',   'h_SIP_lep2',   16, 0, 8.0)
    hst['SIP_max_ll'] = R.TH1D('h_SIP_max_ll', 'h_SIP_max_ll', 16, 0, 8.0)

    hst['pt_log2_sig_MET']    = R.TH1D('h_pt_log2_sig_MET',    'h_pt_log2_sig_MET',    150,  -20,   10)
    hst['pt_diff_MET_pMET']   = R.TH1D('h_pt_diff_MET_pMET',   'h_pt_diff_MET_pMET',   200, -500,  500)
    hst['pt_square_MET_pMET'] = R.TH1D('h_pt_square_MET_pMET', 'h_pt_square_MET_pMET', 200,    0, 1000)

    hst['EF_chH_jetSel']  =  R.TH1D('h_EF_chH_jetSel',   'h_EF_chH_jetSel',   100, 0, 1.0)
    hst['EF_chEm_jetSel'] =  R.TH1D('h_EF_chEm_jetSel',  'h_EF_chEm_jetSel',  100, 0, 1.0)
    hst['EF_chHEm_jetSel'] = R.TH1D('h_EF_chHEm_jetSel', 'h_EF_chHEm_jetSel', 100, 0, 1.0)
    hst['EF_neH_jetSel']  =  R.TH1D('h_EF_neH_jetSel',   'h_EF_neH_jetSel',   100, 0, 1.0)
    hst['EF_neEm_jetSel'] =  R.TH1D('h_EF_neEm_jetSel',  'h_EF_neEm_jetSel',  100, 0, 1.0)
    hst['EF_neHEm_jetSel'] = R.TH1D('h_EF_neHEm_jetSel', 'h_EF_neHEm_jetSel', 100, 0, 1.0)
    hst['EF_chPU_jetSel'] =  R.TH1D('h_EF_chPU_jetSel',  'h_EF_chPU_jetSel',  100, 0, 1.0)
    hst['EF_mu_jetSel']    = R.TH1D('h_EF_mu_jetSel',    'h_EF_mu_jetSel',    100, 0, 1.0)

    hst['pt_jetSel']     = R.TH1D('h_pt_jetSel',     'h_pt_jetSel',     200,    0, 1000)
    hst['pz_jetSel']     = R.TH1D('h_pz_jetSel',     'h_pz_jetSel',    2000, -5000, 5000)
    hst['eta_jetSel']    = R.TH1D('h_eta_jetSel',    'h_eta_jetSel',    200, -5.0,  5.0)
    hst['mass_jetSel']   = R.TH1D('h_mass_jetSel',   'h_mass_jetSel',   200,    0,  100)
    hst['jetId_jetSel']  = R.TH1D('h_jetId_jetSel',  'h_jetId_jetSel',    7, -0.5,  6.5)
    hst['puId_jetSel']   = R.TH1D('h_puId_jetSel',   'h_puId_jetSel',     8, -0.5,  7.5)
    hst['puIdD_jetSel']  = R.TH1D('h_puIdD_jetSel',  'h_puIdD_jetSel',  200,    0,  1.0)
    hst['qgl_jetSel']    = R.TH1D('h_qgl_jetSel',    'h_qgl_jetSel',    100,    0,  1.0)
    hst['flavQG_jetSel'] = R.TH1D('h_flavQG_jetSel', 'h_flavQG_jetSel', 100,    0,  1.0)
    hst['deepB_jetSel']  = R.TH1D('h_deepB_jetSel',  'h_deepB_jetSel',  100,    0,  1.0)
    hst['flavB_jetSel']  = R.TH1D('h_flavB_jetSel',  'h_flavB_jetSel',  100,    0,  1.0)
    hst['nConst_jetSel'] = R.TH1D('h_nConst_jetSel', 'h_nConst_jetSel', 100,    0,  100)

    hst['CoM_log2_pt_lj']    = R.TH1D('h_CoM_log2_pt_lj',    'h_CoM_log2_pt_lj',    150,   -2,  13)
    hst['CoM_log2_pz_lj']    = R.TH1D('h_CoM_log2_pz_lj',    'h_CoM_log2_pz_lj',    150,   -2,  13)
    hst['CoM_log2_p_jet']    = R.TH1D('h_CoM_log2_p_jet',    'h_CoM_log2_p_jet',    150,   -2,  13)
    hst['CoM_dTheta_pz_jet'] = R.TH1D('h_CoM_dTheta_pz_jet', 'h_CoM_dTheta_pz_jet',  64,    0, 3.2)
    hst['CoM_dPhi_pt_jet']   = R.TH1D('h_CoM_dPhi_pt_jet',   'h_CoM_dPhi_pt_jet',    64,    0, 3.2)

    hst['Det_log2_pt_lep']   = R.TH1D('h_Det_log2_pt_lep',   'h_Det_log2_pt_lep',   160,    4,  12)
    hst['Det_log2_pt_jet']   = R.TH1D('h_Det_log2_pt_jet',   'h_Det_log2_pt_jet',   160,    4,  12)
    hst['Det_log2_pt_pMET']  = R.TH1D('h_Det_log2_pt_pMET',  'h_Det_log2_pt_pMET',  160,   -2,   6)
    hst['Det_log2_ptRel_lj'] = R.TH1D('h_Det_log2_ptRel_lj', 'h_Det_log2_ptRel_lj', 250, -2.5, 2.5)
    hst['Det_log2_ptRel_lv'] = R.TH1D('h_Det_log2_ptRel_lv', 'h_Det_log2_ptRel_lv', 200,   -8,   2)

    # hst['mass_min_lep_jet'] = R.TH1D('h_mass_min_lep_jet', 'min mass(lep, jet)', 100, 0, 1000)
    # hst['mass_max_lep_jet'] = R.TH1D('h_mass_max_lep_jet', 'max mass(lep, jet)', 100, 0, 1000)
    # hst['dR_min_lep_jet']   = R.TH1D('h_dR_min_lep_jet',   'min dR(lep, jet)',   100, 0, 10)
    # hst['dR_max_lep_jet']   = R.TH1D('h_dR_max_lep_jet',   'max dR(lep, jet)',   100, 0, 10)
    # hst['dPhi_min_lep_jet'] = R.TH1D('h_dPhi_min_lep_jet', 'min dPhi(lep, jet)',  32, 0, 3.2)
    # hst['dPhi_max_lep_jet'] = R.TH1D('h_dPhi_max_lep_jet', 'max dPhi(lep, jet)',  32, 0, 3.2)
    # hst['dEta_min_lep_jet'] = R.TH1D('h_dEta_min_lep_jet', 'min dEta(lep, jet)', 100, 0, 10)
    # hst['dEta_max_lep_jet'] = R.TH1D('h_dEta_max_lep_jet', 'max dEta(lep, jet)', 100, 0, 10)

    # hst['MT_min_lep_MET'] = R.TH1D('h_MT_min_lep_MET', 'min MT(lep, MET)', 100, 0, 500)
    # hst['MT_max_lep_MET'] = R.TH1D('h_MT_max_lep_MET', 'max MT(lep, MET)', 100, 0, 500)

    # hst['deepB_max_jet'] = R.TH1D('h_deepB_max_jet', 'Maximum AK4 jet deepB tag score', 55, -0.1, 1.0)
    # hst['flavB_max_jet'] = R.TH1D('h_flavB_max_jet', 'Maximum AK4 jet deepFlavB tag score', 55, -0.1, 1.0)
    # hst['nSV_max_jet']   = R.TH1D('h_nSV_max_jet',   'Maximum AK4 jet # of secondary vertices', 7, -0.5, 6.5)

    hst['wgt']      = R.TH1D('h_wgt',      'Event weights (unweighted)', 2000, -10, 10.)
    hst['wgt_wgt']  = R.TH1D('h_wgt_wgt',  'Event weights (weighted)',   2000, -10, 10.)
    hst['nPV']      = R.TH1D('h_nPV',      '# of PVs',      100, -0.5, 99.5)
    hst['nPV_good'] = R.TH1D('h_nPV_good', '# of good PVs', 100, -0.5, 99.5)


    ## 2018 DeepCSV and DeepFlavB cuts: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
    DeepB = {'L': 0.1208, 'M': 0.4168, 'T': 0.7665}
    FlavB = {'L': 0.0490, 'M': 0.2783, 'T': 0.7100}


    ## Loop through events, select, and plot
    nEntries = chains['Events'].GetEntries()
    PS  = ( int(math.floor(1.0*nEntries*nJobs / MAX_EVT)) if MAX_EVT > 0 and MAX_EVT < nEntries*nJobs else 1)
    print '\nEntering loop over %d events, will use prescale factor of %d\n' % (nEntries, PS)

    # print '\nEvent, isSignal, mu_ele, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, mass_llvJ, MT_llvJ, pt_llvJ, eta_llvJ, mass_llJ, pt_jet1, pt_ISRs, dR_max_lep_fat, dPhi_max_lep_fat, mass_max_lep_fat, dR_min_lep_fat, dR_ll_fat, dPhi_ll_fat, dR_lep1_fat, dPhi_lep1_fat, dPhi_llJ_MET, dPhi_llvJ_ISRs, pt_llvJ_ISRs, dPhi_llJ_ISRs, MT_llJ_ISR\n'


    ch = chains['Events']  ## Shortcut expression
    ## Count passing events
    nEvt      = 0
    nPassPre  = 0
    nPassMu   = 0
    nPassMuE  = 0
    nPassLepC = 0
    nPassSel  = 0
    nPassLep   = 0
    # nPassMET   = 0
    # nPassKin   = 0
    # nPassTight = 0
    nPassData  = 0
    # nPassGen   = 0
    for iEvt in range( int(math.ceil(1.0*nEntries / PS)) ):

        if MAX_EVT > 0 and nEvt >= MAX_EVT: break

        ch.GetEntry(iEvt*PS)
        nEvt += 1

        if nEvt % PRT_EVT is 0: print 'Looking at %dst event (#%d / %d)' % (nEvt, iEvt*PS + 1, nEntries)
        # print 'Looking at %dst event (#%d / %d)' % (nEvt, iEvt*PS + 1, nEntries)

        ## Require at least one lepton
        if ch.nMuon + ch.nElectron < 1: continue
        ## If 'twoLep', require at least two leptons
        if 'twoLep' in CAT and ch.nMuon + ch.nElectron < 2: continue

        ## Baseline cuts by category
        if CAT.startswith('SingleMu') and ch.nMuon < 1: continue
        
        ## Remove data events not in Golden JSON
        if isData:
            pass_JSON = EVT_SEL.PassJSON(ch.run, ch.luminosityBlock)
            if not pass_JSON:
                continue
        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            continue

        ## Apply MET cut
        if 'pMET30' in CAT and ch.PuppiMET_pt > 30: continue

        ## Trigger selection in data and MC
        if CAT.startswith('SingleMu') and not (ch.L1_SingleMu22 and (ch.HLT_IsoMu24 or ch.HLT_Mu50)): continue

        nPassPre += 1

        muIdxs  = []  ## Indices of selected muon
        eleIdxs = []  ## Indices of selected electron
        jetIdxs = []  ## Indices of selected AK4 jets

        muVecs   = []  ## TLorentzVectors of selected muons
        muVecTs  = []  ## TLorentzVectors of selected muons
        eleVecs  = []  ## TLorentzVectors of selected electrons
        eleVecTs = []  ## TLorentzVectors of selected electrons
        lepVecs  = []  ## TLorentzVectors of two selected leptons
        lepVecTs = []  ## TLorentzVectors of two selected leptons
        jetVecs  = []  ## TLorentzVectors of selected AK4 jets
        jetVecTs = []  ## TLorentzVectors of selected AK4 jets


        ## Find muon(s) passing cuts
        for iMu in range(ch.nMuon):
            if              ch.Muon_pt [iMu]  < 10.: break
            if          abs(ch.Muon_eta[iMu]) > 2.4: continue
            if   ch.Muon_mediumPromptId[iMu]  != 1 : continue
            ## if     ch.Muon_miniIsoId[iMu]  <  2 : continue  ## Doesn't work!!! - AWB 2021.06.10
            ## Equivalent: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#L699
            if ch.Muon_miniPFRelIso_all[iMu] > 0.20: continue
            if           ch.Muon_mvaTTH[iMu] < -0.4: continue
            if 'lepMVA0p8' in CAT:
                if       ch.Muon_mvaTTH[iMu] <  0.8: continue

            ## Save the selected muon
            muVec  = R.TLorentzVector()
            muVecT = R.TLorentzVector()
            muVec .SetPtEtaPhiM( ch.Muon_pt[iMu], ch.Muon_eta[iMu], ch.Muon_phi[iMu], 0.106 )
            muVecT.SetPtEtaPhiM( ch.Muon_pt[iMu],                0, ch.Muon_phi[iMu], 0.106 )
            muVecs .append(muVec)
            muVecTs.append(muVecT)
            muIdxs .append(iMu)
        ## End loop: for iMu in range(ch.nMuon)

        ## Single muon events must pass HLT_IsoMu24 threshold
        if CAT.startswith('SingleMu') and (len(muVecs) == 0 or muVecs[0].Pt() < 26): continue

        nPassMu += 1
        
        ## Find electron(s) passing cuts
        for iEle in range(ch.nElectron):
            if                 ch.Electron_pt [iEle]  < 10.: break
            if             abs(ch.Electron_eta[iEle]) > 2.5: continue
            ## https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_V2_for_re
            if ch.Electron_mvaFall17V2Iso_WP90[iEle]  != 1 : continue
            if              ch.Electron_mvaTTH[iEle] < -0.4: continue
            if 'lepMVA0p8' in CAT:
                if          ch.Electron_mvaTTH[iEle] <  0.8: continue

            ## Save the selected electron
            eleVec  = R.TLorentzVector()
            eleVecT = R.TLorentzVector()
            eleVec .SetPtEtaPhiM( ch.Electron_pt[iEle], ch.Electron_eta[iEle], ch.Electron_phi[iEle], 0.0005 )
            eleVecT.SetPtEtaPhiM( ch.Electron_pt[iEle],                     0, ch.Electron_phi[iEle], 0.0005 )
            eleVecs .append(eleVec)
            eleVecTs.append(eleVecT)
            eleIdxs .append(iEle)
        ## End loop: for iEle in range(ch.nElectron)

        ## Only keep events with at least one selected lepton
        if len(muVecs) + len(eleVecs) < 1: continue
        ## "oneLep" category requires exactly one pre-selected lepton
        if 'oneLep' in CAT and len(muVecs) + len(eleVecs) > 1: continue
        ## "twoLep" category requires exactly two pre-selected leptons
        if 'twoLep' in CAT and len(muVecs) + len(eleVecs) != 2: continue

        nPassMuE += 1

        ## Require high or low invariant mass lepton pairs
        loM = False
        hiM = False
        for muVec1 in muVecs:
            for muVec2 in muVecs:
                if (muVec1+muVec2).M() > 12 and (muVec1+muVec2).M() < 76: loM = True
                if (muVec1+muVec2).M() > 106: hiM = True
        for eleVec1 in eleVecs:
            for eleVec2 in eleVecs:
                if (eleVec1+eleVec2).M() > 12 and (eleVec1+eleVec2).M() < 76: loM = True
                if (eleVec1+eleVec2).M() > 106: hiM = True
        for muVec in muVecs:
            for eleVec in eleVecs:
                if (muVec+eleVec).M() > 12 and (muVec+eleVec).M() < 76: loM = True
                if (muVec+eleVec).M() > 106: hiM = True
        if 'loM' in CAT and not loM: continue
        if 'hiM' in CAT and not hiM: continue


        ## Define a new lepton category event-by-event
        LepCat = 'NULL'

        if len(eleVecs) == 0:
            if   len(muVecs) == 1: LepCat = 'SingleMu'
            elif len(muVecs) == 2: LepCat = 'DoubleMu'
        elif len(muVecs) == 0:
            if   len(eleVecs) == 1: LepCat = 'SingleEG'
            elif len(eleVecs) == 2: LepCat = 'DoubleEle'
        elif len(muVecs) == 1 and len(eleVecs) == 1:
            if   muVecs[0].Pt() >= eleVecs[0].Pt(): LepCat = 'MuEle'
            else:                                   LepCat = 'EleMu'

        if LepCat == 'NULL':
            continue

        if LepCat == 'SingleMu':
            xMu1 = muIdxs[0]
        if LepCat == 'SingleEG':
            xEle1 = eleIdxs[0]
        if LepCat == 'MuEle' or LepCat == 'EleMu':
            xMu1  = muIdxs[0]
            xEle1 = eleIdxs[0]
        if LepCat == 'DoubleMu':
            xMu1 = muIdxs[0]
            xMu2 = muIdxs[1]
        if LepCat == 'DoubleEle':
            xEle1 = eleIdxs[0]
            xEle2 = eleIdxs[1]

        ## Set lepton vectors by high-pT, low-pT
        if LepCat == 'SingleMu':
            lepVecs  = [muVecs [0]]
            lepVecTs = [muVecTs[0]]
        if LepCat == 'SingleEG':
            lepVecs  = [eleVecs [0]]
            lepVecTs = [eleVecTs[0]]
        if LepCat == 'DoubleMu':
            lepVecs  = [muVecs [0], muVecs [1]]
            lepVecTs = [muVecTs[0], muVecTs[1]]
        if LepCat == 'DoubleEle':
            lepVecs  = [eleVecs [0], eleVecs [1]]
            lepVecTs = [eleVecTs[0], eleVecTs[1]]
        if LepCat == 'MuEle':
            lepVecsTmp  = [muVecs [0], eleVecs [0]]
            lepVecTsTmp = [muVecTs[0], eleVecTs[0]]
        if LepCat == 'EleMu':
            lepVecsTmp  = [eleVecs [0], muVecs [0]]
            lepVecTsTmp = [eleVecTs[0], muVecTs[0]]
        if LepCat != 'SingleMu' and LepCat != 'SingleEle':
            if lepVecs[0].Pt() < lepVecs[1].Pt():
                print '\n*** Super-weird %s event!!! LS = %d, event = %d has pT mis-ordered (%f, %f). ***\n' \
                    % (LepCat, ch.luminosityBlock, ch.event, lepVecs[0].Pt(), lepVecs[1].Pt())


        ## Store MET and PuppiMET vectors (eta depends on leading lepton)
        metVec  = R.TLorentzVector()
        metVecT = R.TLorentzVector()
        metVec .SetPtEtaPhiM( ch.MET_pt, lepVecs[0].Eta(), ch.MET_phi, 0)
        metVecT.SetPtEtaPhiM( ch.MET_pt,                0, ch.MET_phi, 0)

        pMetVec  = R.TLorentzVector()
        pMetVecT = R.TLorentzVector()
        pMetVec .SetPtEtaPhiM( ch.PuppiMET_pt, lepVecs[0].Eta(), ch.PuppiMET_phi, 0)
        pMetVecT.SetPtEtaPhiM( ch.PuppiMET_pt,                0, ch.PuppiMET_phi, 0)

        nPassLepC += 1


        ## Find AK4 jets passing cuts
        for iJet in range(ch.nJet):
            if     ch.Jet_pt [iJet]  <  25: break
            if 'cenJet30' in CAT:
                if ch.Jet_pt [iJet]  <  30: break
            if abs(ch.Jet_eta[iJet]) > 4.7: continue
            if 'cenJet30' in CAT:
                if abs(ch.Jet_eta[iJet]) > 2.4: continue
            ## Tight jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018 (TODO: Update! - AWB 2022.03.11)
            if     ch.Jet_jetId[iJet] <= 1: continue
            ## Loose pileup ID: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID (TODO: Update! - AWB 2022.03.11)
            if      ch.Jet_pt [iJet]  <  50:
                if ch.Jet_puId[iJet]  <=  3: continue
            ## Save 4-vectors of passing AK4 jets
            jetVec  = R.TLorentzVector()
            jetVecT = R.TLorentzVector()
            jetVec .SetPtEtaPhiM( ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
            jetVecT.SetPtEtaPhiM( ch.Jet_pt[iJet],                0, ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
            ## Skip jets which overlap lepton
            lepOverlap = False
            for lepVec in lepVecs:
                if jetVec.DeltaR(lepVec) < 0.4: lepOverlap = True
            if lepOverlap: continue
            ## Save selected AK4 jet
            jetIdxs .append(iJet)
            jetVecs .append(jetVec)
            jetVecTs.append(jetVecT)

        ## Apply jet cuts
        if 'ge1j' in CAT and len(jetVecs) == 0: continue
        if 'eq1j' in CAT and len(jetVecs) != 1: continue
        if 'cenJet30' in CAT:
            jetVecCen = jetVecs[0]
            xCenJ     = jetIdxs[0]

        ## Store secondary vertex 4-vectors
        svVecs = []
        for iSV in range(ch.nSV):
            svVec = R.TLorentzVector()
            svVec.SetPtEtaPhiM( ch.SV_pt[iSV], ch.SV_eta[iSV], ch.SV_phi[iSV], ch.SV_mass[iSV])
            svVecs.append(svVec)


        ############################################
        ###  Start storing kinematic quantities  ###
        ############################################

        ## Shortcut 4-vectors for multi-object quantities
        ## Objs = ['lep1','lep2','MET','lv','lj','ll','jj','ljj','jet1','jet2','jet3','evt']
        lv_vec  = lepVecs [0] + pMetVec 
        lv_vecT = lepVecTs[0] + pMetVecT
        if len(jetVecs) > 0:
            lj_vec  = lepVecs [0] + jetVecs [0] 
            lj_vecT = lepVecTs[0] + jetVecTs[0]
        if len(lepVecs) > 1:
            ll_vec  = lepVecs [0] + lepVecs [1]
            ll_vecT = lepVecTs[0] + lepVecTs[1]
        if len(jetVecs) > 1:
            jj_vec   = jetVecs [0] + jetVecs [1] 
            jj_vecT  = jetVecTs[0] + jetVecTs[1]
            ljj_vec  = lepVecs [0] + jj_vec
            ljj_vecT = lepVecTs[0] + jj_vecT
        ## Create whole event vector
        evt_vec = R.TLorentzVector()
        for lepVec in lepVecs:
            evt_vec = evt_vec + lepVec
        for jetVec in jetVecs:
            evt_vec = evt_vec + jetVec

        nSV_max_jet = 0  ## Maximum number of secondary vertices in one AK4 jet
        for jetVec in jetVecs:
            ## Count matched secondary vertices
            nSV_jet = 0
            for svVec in svVecs:
                if jetVec.DeltaR(svVec) < 0.4: nSV_jet += 1
            nSV_max_jet = max(nSV_max_jet, nSV_jet)


        ##########################
        ## Apply selection cuts ##
        ##########################
        nPassSel += 1
        # ## No medium b-tagged AK4 jet
        # if max_deepB > DeepB['M'] or max_flavB > FlavB['M']: continue
        # nPassBtag += 1
        ## Lepton invariant mass cut to remove low-mass resonances
        if LepCat == 'DoubleMu'  and ll_vec.M() < 12 and ch.Muon_charge    [xMu1] + ch.Muon_charge    [xMu2] == 0: continue
        if LepCat == 'DoubleEle' and ll_vec.M() < 12 and ch.Electron_charge[xMu1] + ch.Electron_charge[xMu2] == 0: continue
        ## For SingleMuon category, muon must be leading lepton
        if CAT.startswith('SingleMu'):
            if not (LepCat == 'SingleMu' or LepCat == 'DoubleMu' or (LepCat == 'MuonEG' and muVecs[0].Pt() > eleVecs[0].Pt())):
                continue
        nPassLep += 1
        # ## HEM veto in data (TODO: Update! - AWB 2022.03.11)
        # if isData and ch.run > 319077 and len(jetVecs) > 0:  ## HEM veto
        #     if jetVecs[0].Eta() < -1.17 and jetVecs[0].Phi() > -1.97 and jetVecs[0].Phi() < -0.47:
        #         continue
        nPassData += 1

        # ## *** Store highly discriminating output variables for MVA algorithm ***
        # ## Event, isSignal, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, pt_llJ, pt_llvJ, MT_llvJ, mass_llJ, mass_llvJ,
        # ## mass_max_lep_fat, mass_min_lep_fat, dR_max_lep_fat, dR_min_lep_fat, dEta_max_lep_fat, dEta_min_lep_fat, dR_ll_fat, dEta_ll_fat'

        # if PRT_MVA_IN:
        #     print '%16d, %d, %d, %6.1f, %6.1f, %6.3f, %6.3f, %6.1f, %6.1f, %6.1f, %5.3f, %6.1f, %6.1f, %6.1f, %5.3f, %5.3f,' \
        #         ' %5.3f, %6.1f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %6.1f, %5.3f, %6.1f' \
        #         % ( ch.luminosityBlock*1000000000 + ch.event, CUT_GEN_BB, (LepCat == 'MuonEG'), ## Event, isSignal, mu_ele
        #             ch.FatJet_msoftdrop[fatIdx], fatVec.M(),                                    ## msoft_fat, mass_fat
        #             max(-0.099, max_deepB), max(-0.099, max_flavB),                             ## deepB_max_jet, flavB_max_jet
        #             llvJ_vec.M(), llvJ_vecT.M(), llvJ_vec.Pt(), abs(llvJ_vec.Eta()),            ## mass_llvJ, MT_llvJ, pt_llvJ, eta_llvJ
        #             llJ_vec.M(), isrVec.Pt(), isrsVec.Pt(),                                     ## mass_llJ, pt_jet1, pt_ISRs
        #             max(lepVecs[0].DeltaR(fatVec), lepVecs[1].DeltaR(fatVec)),                  ## dR_max_lep_fat
        #             min(lepVecs[0].DeltaR(fatVec), lepVecs[1].DeltaR(fatVec)),                  ## dR_min_lep_fat
        #             max(abs(lepVecs[0].DeltaPhi(fatVec)), abs(lepVecs[1].DeltaPhi(fatVec))),    ## dPhi_max_lep_fat
        #             max((lepVecs[0]+fatVec).M(), (lepVecs[1]+fatVec).M()),                      ## mass_max_lep_fat
        #             ll_vec.DeltaR(fatVec), abs(ll_vec.DeltaPhi(fatVec)),                        ## dR_ll_fat, dPhi_ll_fat
        #             lepVecs[0].DeltaR(fatVec), abs(lepVecs[0].DeltaPhi(fatVec)),                ## dR_lep1_fat, dPhi_lep1_fat
        #             abs(llJ_vec.DeltaPhi(metVec)),                                              ## dPhi_llJ_MET
        #             abs(llvJ_vec.DeltaPhi(isrsVec)), (llvJ_vec+isrsVec).Pt(),                   ## dPhi_llvJ_ISRs, pt_llvJ_ISRs
        #             abs(llJ_vec.DeltaPhi(isrsVec)), (llJ_vecT+isrVecT).M() )                    ## dPhi_llJ_ISRs, MT_llJ_ISR
            
            
        #######################
        ## Get event weights ##
        #######################
            
        WGT         = 1.0  ## Overall event weight
        WGT_NO_LUMI = 1.0  ## Event weight before luminosity scaling
        
        if not isData:
            ## Top pT weighting
            if 'TTJets' in in_dir[1]:
                WGT *= EVT_WGT.GetTopPtWgt( top_pt[0], top_pt[1] )
            # ## HEM veto weight (TODO: Update! - AWB 2022.03.11)
            # if len(jetVecs) > 0:
            #     if jetVecs[0].Eta() < -1.17 and jetVecs[0].Phi() > -1.97 and jetVecs[0].Phi() < -0.47:
            #         WGT *= (21090. / LUMI)
            ## GEN negative weight
            if ch.genWeight < 0:
                WGT *= -1
        ## End conditional: if not isData

        WGT_NO_LUMI = WGT  ## Track event weight before cross-section and luminosity scaling
        if not isData:
            WGT *= EVT_WGT.GetXsecPerEvt( in_dir )
            WGT *= LUMI
            ## QCD MuEnriched cross section is 302672 pb, see:
            ## https://indico.cern.ch/event/626095/contributions/2572067/attachments/1452642/2240500/tZq_l_met_170502.pdf
            ## Cross-check also with their W+jets cross section (61526.7) an DYJetsToLL cross sections
        if MAX_EVT > 0 and MAX_EVT < nEntries:
            WGT *= (1.0*nEntries / MAX_EVT)


        #####################
        ## Fill histograms ##
        #####################

        ## Objs = ['lep1','lep2','MET','lv','lj','ll','jj','ljj','jet1','jet2','jet3']
        ## Vars = ['pt','eta','phi','dR','dEta','dPhi','mass','MT','HT']

        ## Figure out which endcap the "whole event" is in
        EC = 1 - 2*(evt_vec.Eta() < 0)

        ## Lepton + jet system, boosted into Center of Mass frame
        lj_boost = lj_vec.BoostVector()
        lep_CoM = R.TLorentzVector()
        jet_CoM = R.TLorentzVector()
        lep_CoM.SetPtEtaPhiM( lepVecs[0].Pt(), lepVecs[0].Eta(), lepVecs[0].Phi(), lepVecs[0].M() )
        jet_CoM.SetPtEtaPhiM( jetVecs[0].Pt(), jetVecs[0].Eta(), jetVecs[0].Phi(), jetVecs[0].M() )
        lep_CoM.Boost(-lj_boost)
        jet_CoM.Boost(-lj_boost)

        ## Fill output tree with variables for NN or BDT training
        if TREE_OUT and len(jetVecs) > 0:
            if LABEL == 'SingleMuon_2018C':    int_vars['proc'][0] =  0
            if LABEL == 'QCD_MuEnriched':      int_vars['proc'][0] =  1
            if LABEL == 'WJetsToLNu':          int_vars['proc'][0] =  2
            if LABEL == 'DYJetsToLL_M-10to50': int_vars['proc'][0] =  3
            if LABEL == 'DYJetsToLL_M-50':     int_vars['proc'][0] =  4
            if LABEL == 'LQsc-SM_ktMin_20':    int_vars['proc'][0] = -1

            flt_vars['wgt'][0] = WGT

            flt_vars['CoM_log2_pt_lj'][0]    = math.log(    lj_vec.Pt(),  2)
            flt_vars['CoM_log2_pz_lj'][0]    = math.log(abs(lj_vec.Pz()), 2)
            flt_vars['CoM_log2_p_jet'][0]    = math.log(jet_CoM.P(), 2)
            flt_vars['CoM_dTheta_pz_jet'][0] = (jet_CoM.Theta() if lj_vec.Pz() > 0 else lep_CoM.Theta())
            flt_vars['CoM_dPhi_pt_jet'][0]   = abs(jet_CoM.DeltaPhi(lj_vec))

            flt_vars['Det_log2_pt_lep'][0]   = math.log(lepVecs[0].Pt(), 2) 
            flt_vars['Det_log2_pt_jet'][0]   = math.log(jetVecs[0].Pt(), 2)
            flt_vars['Det_log2_pt_pMET'][0]  = math.log(   pMetVec.Pt(), 2)
            flt_vars['Det_log2_ptRel_lj'][0] = math.log(jetVecs[0].Pt() / lepVecs[0].Pt(), 2)
            flt_vars['Det_log2_ptRel_lv'][0] = math.log(   pMetVec.Pt() / lepVecs[0].Pt(), 2)

            flt_vars['dPhi_lv'][0]  = abs(lepVecs[0].DeltaPhi(pMetVecT))
            flt_vars['dPhi_lj'][0]  = abs(lepVecs[0].DeltaPhi(jetVecs[0]))
            flt_vars['dPhi_jv'][0]  = abs(jetVecs[0].DeltaPhi(pMetVecT))
            flt_vars['dPhi_ljv'][0] = abs(    lj_vec.DeltaPhi(pMetVecT))

            out_tree.Fill()
            continue  ## Don't fill histograms
        ## End conditional: if TREE_OUT and len(jetVecs) > 0


        ## Whole event variables
        hst['nPV']     .Fill( min( ch.PV_npvs, 100 ), WGT )
        hst['nPV_good'].Fill( min( ch.PV_npvsGood, 100), WGT )
        hst['wgt']     .Fill( max( 0.001, min( 1.999, WGT_NO_LUMI ) ) )
        hst['wgt_wgt'] .Fill( max( 0.001, min( 1.999, WGT_NO_LUMI ) ), WGT_NO_LUMI )

        hst['nLep'].Fill( len(muVecs) + len(eleVecs), WGT)
        hst['nMu'] .Fill( len(muVecs), WGT)
        hst['nEle'].Fill( len(eleVecs), WGT)
        hst['nJet'].Fill( len(jetVecs), WGT)

        hst['pt_evt']  .Fill( evt_vec.Pt(), WGT)
        hst['pz_evt']  .Fill( evt_vec.Pz(), WGT)
        hst['mass_evt'].Fill( evt_vec.M(),  WGT)

        ## pt and pz of lepton + jet system in the detector rest frame
        hst['CoM_log2_pt_lj'].Fill( math.log(    lj_vec.Pt(),  2), WGT)
        hst['CoM_log2_pz_lj'].Fill( math.log(abs(lj_vec.Pz()), 2), WGT)
        ## Momentum of jet and lepton in CoM frame are the same, roughly half the mass
        hst['CoM_log2_p_jet'].Fill( math.log(jet_CoM.P(), 2), WGT)
        ## Use theta of jet relative to original 'z' direction of lepton + jet system
        hst['CoM_dTheta_pz_jet'].Fill( (jet_CoM.Theta() if lj_vec.Pz() > 0 else lep_CoM.Theta()), WGT)
        ## Use phi of jet relative to original phi direction of lepton + jet system
        hst['CoM_dPhi_pt_jet'].Fill( abs(jet_CoM.DeltaPhi(lj_vec)), WGT)
        ## Also fill some detector-frame quantities with logs
        hst['Det_log2_pt_lep']  .Fill( math.log(lepVecs[0].Pt(), 2), WGT) 
        hst['Det_log2_pt_jet']  .Fill( math.log(jetVecs[0].Pt(), 2), WGT)
        hst['Det_log2_pt_pMET'] .Fill( math.log(   pMetVec.Pt(), 2), WGT)
        hst['Det_log2_ptRel_lj'].Fill( math.log(jetVecs[0].Pt() / lepVecs[0].Pt(), 2), WGT)
        hst['Det_log2_ptRel_lv'].Fill( math.log(   pMetVec.Pt() / lepVecs[0].Pt(), 2), WGT)
        

        ## Lepton properties
        if LepCat == 'SingleMu' or LepCat == 'DoubleMu' or (LepCat == 'MuonEG' and muVecs[0].Pt() > eleVecs[0].Pt()):
            hst['charge_lep1'] .Fill( ch.Muon_charge[xMu1], WGT)
            hst['flavor_lep1'] .Fill( 1, WGT)
            hst['mvaTTH_lep1'] .Fill( ch.Muon_mvaTTH[xMu1], WGT)
            hst['miniIso_lep1'].Fill( ch.Muon_miniPFRelIso_all[xMu1], WGT)
            hst['SIP_lep1']    .Fill( min( abs(ch.Muon_sip3d[xMu1]), 7.99), WGT)
            hst['ID_lep1']     .Fill(   ch.Muon_mediumPromptId[xMu1] + \
                                      2*ch.Muon_tightId       [xMu1], WGT)

            if LepCat == 'SingleMu':
                hst['charge_sum_ll'].Fill( ch.Muon_charge[xMu1], WGT)
                hst['flavor_sum_ll'].Fill( 1, WGT)

            if LepCat == 'DoubleMu':
                hst['charge_lep2']  .Fill( ch.Muon_charge[xMu2], WGT)
                hst['charge_sum_ll'].Fill( ch.Muon_charge[xMu1] + ch.Muon_charge[xMu2], WGT)

                hst['flavor_lep2']  .Fill( 1, WGT)
                hst['flavor_sum_ll'].Fill( 2, WGT)

                hst['mvaTTH_lep2']  .Fill( ch.Muon_mvaTTH[xMu2], WGT)
                hst['mvaTTH_min_ll'].Fill( min(ch.Muon_mvaTTH[xMu1], ch.Muon_mvaTTH[xMu2]), WGT)

                hst['miniIso_lep2']  .Fill( ch.Muon_miniPFRelIso_all[xMu2], WGT)
                hst['miniIso_max_ll'].Fill( min( max(ch.Muon_miniPFRelIso_all[xMu1], ch.Muon_miniPFRelIso_all[xMu2]), 0.79), WGT)

                hst['SIP_lep2']  .Fill( min( abs(ch.Muon_sip3d[xMu2]), 7.99), WGT)
                hst['SIP_max_ll'].Fill( min( max(abs(ch.Muon_sip3d[xMu1]), abs(ch.Muon_sip3d[xMu2])), 7.99), WGT)

                hst['ID_lep2'] .Fill(   ch.Muon_mediumPromptId[xMu2] + \
                                      2*ch.Muon_tightId       [xMu2], WGT)
            ## End conditional: if LepCat == 'DoubleMu'

            if LepCat == 'MuonEG' and muVecs[0].Pt() > eleVecs[0].Pt():
                hst['charge_lep2']  .Fill( ch.Electron_charge[xEle1], WGT)
                hst['charge_sum_ll'].Fill( ch.Muon_charge[xMu1] + ch.Electron_charge[xEle1], WGT)

                hst['flavor_lep2']  .Fill( -1, WGT)
                hst['flavor_sum_ll'].Fill(  0, WGT)

                hst['mvaTTH_lep2']  .Fill( ch.Electron_mvaTTH[xEle1], WGT)
                hst['mvaTTH_min_ll'].Fill( min(ch.Muon_mvaTTH[xMu1], ch.Electron_mvaTTH[xEle1]), WGT)

                hst['miniIso_lep2']  .Fill( ch.Electron_miniPFRelIso_all[xEle1], WGT)
                hst['miniIso_max_ll'].Fill( min( max(ch.Muon_miniPFRelIso_all[xMu1], ch.Electron_miniPFRelIso_all[xEle1]), 0.79), WGT)

                hst['SIP_lep2']  .Fill( min( abs(ch.Electron_sip3d[xEle1]), 7.99), WGT)
                hst['SIP_max_ll'].Fill( min( max(abs(ch.Muon_sip3d[xMu1]), abs(ch.Electron_sip3d[xEle1])), 7.99), WGT)

                hst['ID_lep2'].Fill(   ch.Electron_mvaFall17V2Iso_WPL   [xEle1] + \
                                     2*ch.Electron_mvaFall17V2noIso_WP90[xEle1] + \
                                     4*ch.Electron_mvaFall17V2Iso_WP90  [xEle1] + \
                                     8*ch.Electron_mvaFall17V2noIso_WP80[xEle1] + \
                                    16*ch.Electron_mvaFall17V2Iso_WP80  [xEle1], WGT)
            ## End conditional: if LepCat == 'MuonEG' and muVecs[0].Pt() > eleVecs[0].Pt()

        ## End conditional: if LepCat == 'SingleMu' or LepCat == 'DoubleMu' or (LepCat == 'MuonEG' and muVecs[0].Pt() > eleVecs[0].Pt()):

        ## Objs = ['lep1','lep2','MET','lv','lj','ll','jj','ljj','jet1','jet2','jet3']
        ## Vars = ['pt','eta','phi','dR','dEta','dPhi','mass','MT','HT']

        ## Leading lepton
        hst['pt_lep1'] .Fill( lepVecs[0].Pt(),     WGT)
        hst['pz_lep1'] .Fill( lepVecs[0].Pz()*EC,  WGT)
        hst['eta_lep1'].Fill( lepVecs[0].Eta()*EC, WGT)
        hst['phi_lep1'].Fill( lepVecs[0].Phi(),    WGT)

        ## MET
        hst['pt_MET']  .Fill( metVec.Pt(),   WGT)
        hst['phi_MET'] .Fill( metVec.Phi(),  WGT)
        hst['pt_pMET'] .Fill( pMetVec.Pt(),  WGT)
        hst['phi_pMET'].Fill( pMetVec.Phi(), WGT)

        hst['pt_log2_sig_MET']   .Fill( math.log(ch.MET_significance, 2), WGT)
        hst['pt_diff_MET_pMET']  .Fill( metVec.Pt() - pMetVec.Pt(),     WGT)
        hst['pt_square_MET_pMET'].Fill( math.sqrt( pow(metVec.Pt(), 2) + pow(pMetVec.Pt(), 2)), WGT)

        ## Jets
        for i in range(min(len(jetVecs),3)):
            hst['pt_jet%d'  % (i+1)].Fill( jetVecs[i].Pt(),     WGT)
            hst['pz_jet%d'  % (i+1)].Fill( jetVecs[i].Pz()*EC,  WGT)
            hst['eta_jet%d' % (i+1)].Fill( jetVecs[i].Eta()*EC, WGT)
            hst['phi_jet%d' % (i+1)].Fill( jetVecs[i].Phi(),    WGT)

        ## Leading jet
        if len(jetVecs) > 0:
            xJ = (xCenJ if 'cenJet' in CAT else jetIdxs[0])
            jetVecSel = (jetVecCen if 'cenJet' in CAT else jetVecs[0])

            hst['EF_chH_jetSel']  .Fill( ch.Jet_chHEF   [xJ], WGT)
            hst['EF_chEm_jetSel'] .Fill( ch.Jet_chEmEF  [xJ], WGT)
            hst['EF_chHEm_jetSel'].Fill( ch.Jet_chHEF   [xJ] + ch.Jet_chEmEF[xJ], WGT)
            hst['EF_neH_jetSel']  .Fill( ch.Jet_neHEF   [xJ], WGT)
            hst['EF_neEm_jetSel'] .Fill( ch.Jet_neEmEF  [xJ], WGT)
            hst['EF_neHEm_jetSel'].Fill( ch.Jet_neHEF   [xJ] + ch.Jet_neEmEF[xJ], WGT)
            hst['EF_chPU_jetSel'] .Fill( ch.Jet_chFPV0EF[xJ], WGT)
            hst['EF_mu_jetSel']   .Fill( ch.Jet_muEF    [xJ], WGT)
            
            hst['pt_jetSel']    .Fill( ch.Jet_pt            [xJ], WGT)
            hst['pz_jetSel']    .Fill( jetVecSel.Pz()*EC,         WGT)
            hst['eta_jetSel']   .Fill( ch.Jet_eta[xJ]*EC,         WGT)
            hst['mass_jetSel']  .Fill( ch.Jet_mass          [xJ], WGT)
            hst['jetId_jetSel'] .Fill( ch.Jet_jetId         [xJ], WGT)
            hst['puId_jetSel']  .Fill( ch.Jet_puId          [xJ], WGT)
            hst['puIdD_jetSel'] .Fill( ch.Jet_puIdDisc      [xJ], WGT)
            hst['qgl_jetSel']   .Fill( ch.Jet_qgl           [xJ], WGT)
            hst['flavQG_jetSel'].Fill( ch.Jet_btagDeepFlavQG[xJ], WGT)
            hst['deepB_jetSel'] .Fill( ch.Jet_btagDeepB     [xJ], WGT)
            hst['flavB_jetSel'] .Fill( ch.Jet_btagDeepFlavB [xJ], WGT)
            hst['nConst_jetSel'].Fill( ch.Jet_nConstituents [xJ], WGT)

        ## Leading lepton and MET
        hst['MT_lv']  .Fill( lv_vecT.M(),                            WGT)
        hst['HT_lv']  .Fill(     lepVecTs[0].Pt() +   pMetVecT.Pt(), WGT)
        hst['dPhi_lv'].Fill( abs(lepVecTs[0].DeltaPhi(pMetVecT)),    WGT)

        ## Leading lepton and leading jet
        if len(jetVecs) > 0:
            hst['pt_lj']  .Fill( lj_vec.Pt(),     WGT)
            hst['pz_lj']  .Fill( lj_vec.Pz()*EC,  WGT)
            hst['eta_lj'] .Fill( lj_vec.Eta()*EC, WGT)
            hst['mass_lj'].Fill( lj_vec.M(),      WGT)

            hst['HT_lj']   .Fill(     lepVecs[0].Pt() +   jetVecs[0].Pt(),   WGT)
            hst['dR_lj']   .Fill(     lepVecs[0].DeltaR  (jetVecs[0]),       WGT)
            hst['dPhi_lj'] .Fill( abs(lepVecs[0].DeltaPhi(jetVecs[0])),      WGT)
            hst['dEta_lj'] .Fill( abs(lepVecs[0].Eta() -  jetVecs[0].Eta()), WGT)
            hst['dPhi_jv'] .Fill( abs(jetVecs[0].DeltaPhi(pMetVecT)),        WGT)
            hst['dPhi_ljv'].Fill( abs(    lj_vec.DeltaPhi(pMetVecT)),        WGT)

        ## Trailing lepton and ll pair
        if len(lepVecs) > 1:
            hst['pt_lep2'] .Fill( lepVecs[1].Pt(),     WGT)
            hst['pz_lep2'] .Fill( lepVecs[1].Pz()*EC,  WGT)
            hst['eta_lep2'].Fill( lepVecs[1].Eta()*EC, WGT)
            hst['phi_lep2'].Fill( lepVecs[1].Phi(),    WGT)

            hst['pt_ll']  .Fill(     ll_vec.Pt(),        WGT)
            hst['pz_ll']  .Fill(     ll_vec.Pz()*EC,     WGT)
            hst['eta_ll'] .Fill(     ll_vec.Eta()*EC,    WGT)
            hst['mass_ll'].Fill( min(ll_vec.M(), 399.9), WGT)

            hst['HT_ll']  .Fill(     lepVecs[0].Pt() +   lepVecs[1].Pt(),   WGT)
            hst['dR_ll']  .Fill(     lepVecs[0].DeltaR  (lepVecs[1]),       WGT)
            hst['dPhi_ll'].Fill( abs(lepVecs[0].DeltaPhi(lepVecs[1])),      WGT)
            hst['dEta_ll'].Fill( abs(lepVecs[0].Eta() -  lepVecs[1].Eta()), WGT)


    ## End loop: for iEvt in range(nEntries):

    print '\nFinished loop over events'
    
    print '\nOut of %d events, %d pass pre-selection, %d muon (%d +ele) cuts, %d LepCat, %d selection, %d lepton cuts, %d data' \
        % (nEvt, nPassPre, nPassMu, nPassMuE, nPassLepC, nPassSel, nPassLep, nPassData)

    out_file.cd()

    for key in hst.keys():
        if hst[key].Integral() == 0:
            del hst[key]
            continue
        hst[key].SetLineColor(R.kBlack)
        hst[key].SetLineWidth(2)
        
    print '\nSaved histograms to output file %s\n' % out_file_str

    out_file.Write('', R.TObject.kOverwrite)
    out_file.Close()


if __name__ == '__main__':
    main()
