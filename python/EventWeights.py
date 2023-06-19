

################################################################################
###  Functions to weight events according to luminosity, pileup, and LHE HT  ###
################################################################################

import math

import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
VERBOSE = False

class GetWeights:

    def __init__(self):
        print '\nInside "init" for EventWeights (does nothing right now)'
        

    ## Reweight ttbar events based on top pT
    def GetTopPtWgt(self, top_pt, tbar_pt):
        ## From https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat
        ## Not clear whether this is actually appropriate for MadgraphMLM ...
        return math.sqrt( math.exp(0.0615 - 0.0005*top_pt) * math.exp(0.0615 - 0.0005*tbar_pt) )
    ## End function: GetTopPtWgt(self, top_pt, tbar_pt)


    ## Compute Xsec/nEvt ratio for each sample
    def GetXsecPerEvt(self, in_dir):
        ## Cross sections (in pb) currently from:
        ## https://indico.cern.ch/event/626095/contributions/2572067/attachments/1452642/2240500/tZq_l_met_170502.pdf
        ## Event numbers from DAS: dataset=/[SAMPLE_NAME]/*UL18NanoAODv9*/NANOAODSIM

        if 'QCD/QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8'          in in_dir: return ( 302672.0   / 17392472.)
        if 'WJets/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'          in in_dir: return (  61526.7   / 81051269.)
        if 'ZJets/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8' in in_dir: return (  22635.1   / 94452816.)
        if 'ZJets/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'     in in_dir: return (   6025.2   / 96233328.)
        if 'LeptonJet/LQsc-SM_ktMin_20'                                  in in_dir: return (     26.61  /   149000.)
        if 'LeptonJet/LQsc-SM_up_mu_ktMin_20'                            in in_dir: return (     26.61  /     1000.)
        if 'LeptonJet/LQsc-SM_up_mu_ktMin_60'                            in in_dir: return (      1.333 /     1000.)

        print 'Really ought to quit right away, trying to return sqrt(-1) ...\n\n'
        return math.sqrt(-1)
    ## End function: GetXsecPerEvt(self, in_dir)


## End class GetWeights


