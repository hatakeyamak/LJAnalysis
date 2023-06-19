#! /usr/bin/env python

## First run 'conda activate mlenv'

## Plots for lepton-quark scattering to lepton + jet + no MET final state

## Coffea implementation of Lepton + Jet analysis code (LeptonJet.py), following:
## https://indico.cern.ch/event/1129899/contributions/4741912/attachments/2397529/4099651/CoffeaATTF_LindseyGray_25022022.pdf
## https://github.com/nsmith-/TTGamma_LongExercise/blob/FullAnalysis/ttgamma/processor.py

import awkward as ak
import numpy as np
import hist
from hist import Hist
import matplotlib.pyplot as plt
from coffea import nanoevents, processor
## from coffea import hist as coffea_hist
from coffea.lumi_tools import LumiMask
from coffea.nanoevents.methods import vector
#KH import uproot3
import uproot

import os
import sys
import pickle
import subprocess

LABEL   = 'XXX'
CAT     = 'SingleMu_oneLep_lepMVA0p8_pMET30_ge1j'
# CAT     = 'SingleMu_oneLep_lepMVA0p8_pMET30_eq1j'
JSON    = 'data/JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
NJOB    = 1
IJOB    = 0
TEST    = True
VERBOSE = False

class LeptonJetProcessor(processor.ProcessorABC):
    ## Plot properties of lepton + jet events

    ## Initialize options
    def __init__(self):
        print('Inside LeptonJetProcessor __init__')

        self.cat    = CAT
        self.json   = LumiMask(JSON)
        self.isData = ('SingleMuon' in LABEL or 'SingleElectron' in LABEL)
        self.verb   = VERBOSE
        if self.verb:
            print('Set cat = %s, isData = %d, got json from:' % (self.cat, self.isData))
            print(JSON)
        
        self.mu_pt      =  10.
        self.ele_pt     =  10.
        self.mu_trg_pt  =  26.
        self.ele_trg_pt =  35.
        self.mu_eta     =  2.4
        self.ele_eta    =  2.5
        self.mu_iso     =  0.2
        self.mu_MVA     = -0.4 if not 'lepMVA0p8' in self.cat else 0.8
        self.ele_MVA    = -0.4 if not 'lepMVA0p8' in self.cat else 0.8
        self.jet_pt     =  25. if not 'cenJet30'  in self.cat else 30.
        self.jet_eta    =  4.7 if not 'cenJet30'  in self.cat else 2.4
    ## End __init__


    ## Main data-processing and histogram-making function
    def process(self, evt):
        print('Inside LeptonJetProcessor process, total of %d events' % len(evt))

        ## print(evt.fields) to print all properties of 'events'

        ## Require at least one lepton
        evt = evt[ak.num(evt.Muon) + ak.num(evt.Electron) > 0]
        if self.verb:
            print('After >= 1 lepton cut, have %d events' % len(evt))
        ## If 'twoLep', require at least two leptons
        if 'twoLep' in self.cat:
            evt = evt[ak.num(evt.Muon) + ak.num(evt.Electron) > 0]
            if self.verb:
                print('After >= 2 lepton cut, have %d events' % len(evt))
        
        ## Baseline cuts by category
        if self.cat.startswith('SingleMu'):
            evt = evt[ak.num(evt.Muon) > 0]
            if self.verb:
                print('After >= 1 muon cut, have %d events' % len(evt))
            evt = evt[evt.L1.SingleMu22 & (evt.HLT.IsoMu24 | evt.HLT.Mu50)]
            if self.verb:
                print('After trigger cut, have %d events' % len(evt))

        ## Remove data evt not in Golden JSON
        if self.isData:
            evt = evt[self.json(evt.run, evt.luminosityBlock)]
            if self.verb:
                print('After JSON cut, have %d events' % len(evt))

        ## Remove events with no good collision vertices
        evt = evt[evt.PV.npvsGood > 0]
        if self.verb:
            print('After nPV cut, have %d events' % len(evt))

        ## Apply MET cut
        if 'pMET30' in self.cat:
            evt = evt[evt.PuppiMET.pt <= 30]
            if self.verb:
                print('After MET < 30 cut, have %d events' % len(evt))

        
        ############################
        ## Find muon(s) passing cuts
        ############################
        evt['selMu'] = evt.Muon[ (
            (    evt.Muon.pt               >= self.mu_pt  ) &
            (abs(evt.Muon.eta)             <= self.mu_eta ) &
            (    evt.Muon.mediumPromptId   == 1           ) &
            (    evt.Muon.miniPFRelIso_all <= self.mu_iso ) &
            (    evt.Muon.mvaTTH           >= self.mu_MVA ) ) ]

        ## Require at least one selected muon
        if self.cat.startswith('SingleMu'):
            evt = evt[ak.num(evt.selMu) > 0]
            if self.verb:
                print('After SingleMu cut, have %d events' % len(evt))
            ## Apply further cuts on the leading muon
            evt = evt[ak.firsts(evt.selMu).pt >= self.mu_trg_pt]
            if self.verb:
                print('After mu_trg_pt cut, have %d events' % len(evt))


        ################################
        ## Find electron(s) passing cuts
        ################################
        evt['selEle'] = evt.Electron[ (
            (    evt.Electron.pt                  >= self.ele_pt  ) &
            (abs(evt.Electron.eta)                <= self.ele_eta ) &
            (    evt.Electron.mvaFall17V2Iso_WP90 == 1            ) &
            (    evt.Electron.mvaTTH              >= self.ele_MVA ) ) ]

        ## Require at least one selected electron
        if self.cat.startswith('SingleEle'):
            evt = evt[ak.num(evt.selEle) > 0]
            if self.verb:
                print('After SingleEle cut, have %d events' % len(evt))
            ## Apply further cuts on the leading electron
            evt = evt[ak.firsts(evt.selEle).pt >= self.ele_trg_pt]
            if self.verb:
                print('After ele_trg_pt cut, have %d events' % len(evt))
        

        ## Only keep evt with at least one selected lepton
        evt = evt[ak.num(evt.selMu) + ak.num(evt.selEle) > 0]
        if self.verb:
            print('After >= 1 selLep cut, have %d events' % len(evt))
        ## "oneLep" and "twoLep" categories requires exactly 1 or 2 pre-selected leptons
        if 'oneLep' in self.cat:
            evt = evt[ak.num(evt.selMu) + ak.num(evt.selEle) == 1]
            if self.verb:
                print('After oneLep cut, have %d events' % len(evt))
        if 'twoLep' in self.cat:
            evt = evt[ak.num(evt.selMu) + ak.num(evt.selEle) == 2]
            if self.verb:
                print('After twoLep cut, have %d events' % len(evt))

        ## Set muon and electron masses
        evt['selMu',  'mass'] = 0.10566
        evt['selEle', 'mass'] = 0.000511
        ## Set muon and electron 4-vectors
        ## https://coffeateam.github.io/coffea/modules/coffea.nanoevents.methods.vector.html
        evt['selMu', 'p4'] = ak.zip( { 'pt':   evt.selMu.pt,
                                       'eta':  evt.selMu.eta,
                                       'phi':  evt.selMu.phi,
                                       'mass': evt.selMu.mass },
                                     with_name='PtEtaPhiMLorentzVector',
                                     behavior=vector.behavior )
        evt['selEle', 'p4'] = ak.zip( { 'pt':   evt.selEle.pt,
                                        'eta':  evt.selEle.eta,
                                        'phi':  evt.selEle.phi,
                                        'mass': evt.selEle.mass },
                                      with_name='PtEtaPhiMLorentzVector',
                                      behavior=vector.behavior )
                                     
        ## Combine muons and electrons into leptons
        evt['selLep'] = ak.concatenate([evt.selMu, evt.selEle], axis=1)
        ## Construct di-lepton pairs
        if 'twoLep' in self.cat:
            evt['pair']   = ak.argcombinations(evt.selLep, 2, fields=['l1','l2'])
            evt['pairLo'] = evt.pair[ ( ((evt.pair.l1 + evt.pair.l2).mass >=  12) &
                                        ((evt.pair.l1 + evt.pair.l2).mass <=  76) ) ]
            evt['pairHi'] = evt.pair[ (  (evt.pair.l1 + evt.pair.l2).mass >= 106  ) ]
            if 'loM' in self.cat:
                evt = evt[ak.num(evt.pairLo) > 0]
            if 'hiM' in self.cat:
                evt = evt[ak.num(evt.pairHi) > 0]

        # ## Define a new lepton category event-by-event
        # evt['LepCat'] = 'NULL'
        # evt['LepCat'][( (ak.num(evt.selMu) == 1) & (ak.num(evt.selEle) == 0) )] = 'SingleMu' 
        # evt['LepCat'][( (ak.num(evt.selMu) == 0) & (ak.num(evt.selEle) == 1) )] = 'SingleEle' 
        # evt['LepCat'][( (ak.num(evt.selMu) == 2) & (ak.num(evt.selEle) == 0) )] = 'DoubleMu' 
        # evt['LepCat'][( (ak.num(evt.selMu) == 0) & (ak.num(evt.selEle) == 2) )] = 'DoubleEle' 
        # evt['LepCat'][( (ak.num(evt.selMu) == 1) & (ak.num(evt.selEle) == 1) &
        #                 (ak.firsts(evt.selMu).pt >= ak.firsts(evt.selEle).pt) )] = 'MuEle'
        # evt['LepCat'][( (ak.num(evt.selMu) == 1) & (ak.num(evt.selEle) == 1) &
        #                 (ak.firsts(evt.selMu).pt <  ak.firsts(evt.selEle).pt) )] = 'EleMu'
        # evt = evt[evt.LepCat != 'NULL']


        ###########################
        ## Find jet(s) passing cuts
        ###########################
        if self.verb:
            print('Before jet logic, have %d events' % len(evt))

        evt['selJet'] = evt.Jet[ (
            (    evt.Jet.pt    >= self.jet_pt         ) &
            (abs(evt.Jet.eta)  <= self.jet_eta        ) &
            (    evt.Jet.jetId >= 2                   ) &
            (    evt.Jet.puId  >= 4*(evt.Jet.pt < 50) ) ) ]
        ## Require at least one selected jet
        evt = evt[ak.num(evt.selJet) > 0]
        if self.verb:
            print('After >= 1 jet cut, have %d events' % len(evt))
        ## Require jets to not overlap any selected leptons
        evt['selJet'] = evt.selJet[ evt.selJet.delta_r(evt.selJet.nearest(evt.selLep)) >= 0.4 ]

        ## Set jet 4-vectors
        evt['selJet', 'p4'] = ak.zip( { 'pt':   evt.selJet.pt,
                                        'eta':  evt.selJet.eta,
                                        'phi':  evt.selJet.phi,
                                        'mass': evt.selJet.mass },
                                      with_name='PtEtaPhiMLorentzVector',
                                      behavior=vector.behavior )

        ## Apply event-level jet cuts
        if 'ge1j' in self.cat:
            evt = evt[ak.num(evt.selJet) > 0]
            if self.verb:
                print('After ge1j cut, have %d events' % len(evt))
        if 'eq1j' in self.cat:
            evt = evt[ak.num(evt.selJet) == 1]
            if self.verb:
                print('After eq1j cut, have %d events' % len(evt))

        if self.verb:
            print('Before filling histograms, have %d events' % len(evt))

        ## https://hist.readthedocs.io/en/latest/user-guide/quickstart.html
        hists = []
        ## Whole event variables
        hists.append( Hist.new.Reg( 100, -0.5, 99.5, name='nPV'     ).Double().fill( evt.PV.npvs ) )
        hists.append( Hist.new.Reg( 100, -0.5, 99.5, name='nPV_good').Double().fill( evt.PV.npvsGood ) )
        hists.append( Hist.new.Reg( 2000, -10, 10.,  name='wgt'     ).Double().fill( np.zeros_like(evt.PV.npvs) + 1 ) )
        hists.append( Hist.new.Reg( 2000, -10, 10.,  name='wgt_wgt' ).Double().fill( np.zeros_like(evt.PV.npvs) + 1 ) )

        hists.append( Hist.new.Reg(  5, -0.5, 4.5, name='nLep').Double().fill( ak.num(evt.selLep) ) )
        hists.append( Hist.new.Reg(  5, -0.5, 4.5, name='nMu' ).Double().fill( ak.num(evt.selMu)  ) )
        hists.append( Hist.new.Reg(  5, -0.5, 4.5, name='nEle').Double().fill( ak.num(evt.selEle) ) )
        hists.append( Hist.new.Reg( 10, -0.5, 9.5, name='nJet').Double().fill( ak.num(evt.selJet) ) )

        full_evt_vec = vector.LorentzVector.sum(evt.selLep.p4) + vector.LorentzVector.sum(evt.selJet.p4)
        hists.append( Hist.new.Reg(  200,     0, 1000, name='pt_evt'  ).Double().fill( full_evt_vec.pt ) )
        hists.append( Hist.new.Reg( 2000, -5000, 5000, name='pz_evt'  ).Double().fill( full_evt_vec.pz ) )
        hists.append( Hist.new.Reg(  200,     0, 2000, name='mass_evt').Double().fill( full_evt_vec.mass ) )

        lep_vec = ak.firsts(evt.selLep.p4)
        jet_vec = ak.firsts(evt.selJet.p4)
        lj_vec  = lep_vec + jet_vec
        lep_CoM = lep_vec.boost(-lj_vec.boostvec)
        jet_CoM = jet_vec.boost(-lj_vec.boostvec)
        ## pt and pz of lepton + jet system in the detector rest frame
        hists.append( Hist.new.Reg( 150,   -2,  13, name='CoM_log2_pt_lj').Double().fill( np.log2(lj_vec.pt) ) )
        hists.append( Hist.new.Reg( 150,   -2,  13, name='CoM_log2_pz_lj').Double().fill( np.log2(np.abs(lj_vec.pz)) ) )
        ## Momentum of jet and lepton in CoM frame are the same, roughly half the mass
        hists.append( Hist.new.Reg( 150,   -2,  13, name='CoM_log2_p_jet').Double().fill( np.log2(jet_CoM.p) ) )
        ## Use theta of jet relative to original 'z' direction of lepton + jet system
        hists.append( Hist.new.Reg(  64,    0, 3.2, name='CoM_dTheta_pz_jet').Double().fill(        jet_CoM.theta*(lj_vec.pz > 0) + lep_CoM.theta*(lj_vec.pz <= 0)  ) )
        hists.append( Hist.new.Reg( 100,   -1,   1, name='CoM_cosTheta_jet' ).Double().fill( np.cos(jet_CoM.theta*(lj_vec.pz > 0) + lep_CoM.theta*(lj_vec.pz <= 0)) ) )
        ## Use phi of jet relative to original phi direction of lepton + jet system
        hists.append( Hist.new.Reg(  64,    0, 3.2, name='CoM_dPhi_pt_jet').Double().fill(        np.abs(jet_CoM.delta_phi(lj_vec))  ) )
        hists.append( Hist.new.Reg( 100,   -1,   1, name='CoM_cosPhi_jet' ).Double().fill( np.cos(np.abs(jet_CoM.delta_phi(lj_vec))) ) )
        ## Also fill some detector-frame quantities with logs
        hists.append( Hist.new.Reg( 160,    4,  12, name='Det_log2_pt_lep'  ).Double().fill( np.log2(lep_vec.pt) ) )
        hists.append( Hist.new.Reg( 160,    4,  12, name='Det_log2_pt_jet'  ).Double().fill( np.log2(jet_vec.pt) ) )
        hists.append( Hist.new.Reg( 160,   -2,   6, name='Det_log2_pt_pMET' ).Double().fill( np.log2(evt.PuppiMET.pt) ) )
        hists.append( Hist.new.Reg( 250, -2.5, 2.5, name='Det_log2_ptRel_lj').Double().fill( np.log2(     jet_vec.pt / lep_vec.pt) ) )
        hists.append( Hist.new.Reg( 200,   -8,   2, name='Det_log2_ptRel_lv').Double().fill( np.log2(evt.PuppiMET.pt / lep_vec.pt) ) )

        # ## hists.append( Hist.new.Reg( , , , name='').Double().fill(  ) )
        

        return hists


    
    def postprocess(self, accumulator):
        print('Inside LeptonJetProcessor postprocess')
        return accumulator


if __name__ == "__main__":

    print('\nInside LeptonJetC.py, about to define input options')

    if len(sys.argv) > 1:
        print('\nLABEL changed from %s to %s' % (LABEL, str(sys.argv[1])))
        LABEL = str(sys.argv[1])
    if len(sys.argv) > 2:
        print('\nCAT changed from %s to %s' % (CAT, str(sys.argv[2])))
        CAT = str(sys.argv[2])
    if len(sys.argv) > 3:
        print('\nMAX_EVT changed from %d to %d' % (MAX_EVT, int(sys.argv[3])))
        MAX_EVT = int(sys.argv[3])
        PRT_EVT = MAX_EVT / 100
    if len(sys.argv) > 4:
        print('\nWill split into %d jobs, run job #%d' % (int(sys.argv[4]), int(sys.argv[5])))
        NJOB = int(sys.argv[4])
        IJOB = int(sys.argv[5])

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
        print('\n\n***  TRIED TO APPLY LABEL %s TO DIRECTORY %s  ***' % (LABEL, in_dir))
        print('***  CORRECT? QUITTING!  ***')
        sys.exit()

    nFiles = len(subprocess.check_output(['ls', in_dir]).splitlines()) - 1
    nJobs  = min(NJOB, nFiles)
    print('\n\n*** Found %d total input files ***\n\n' % nFiles)
    if IJOB >= nJobs:
        print('\n\n***  TRIED TO RUN IJOB = %d WITH ONLY %d JOBS TOTAL!!! ***' % (IJOB, nJobs))
        sys.exit()
    iFile = 0
    for f_name in subprocess.check_output(['ls', in_dir]).splitlines():
        if TEST: continue
        f_name = str(f_name).replace("b'","").replace("root'","root")
        if not '.root' in str(f_name): continue
        if 'LQsc-SM_up_mu' in str(f_name) and not 'GTfix' in str(f_name): continue
        if (iFile % nJobs) == IJOB:
            in_file_names.append(in_dir+str(f_name))
            print('Appending file: %s' % in_file_names[-1])
        iFile += 1
    if TEST:
        ## Large (3M evt), medium (350k), and small (20k) files for testing
        in_file_names = ['/cms/data/store/user/abrinke1/NanoAOD/2018/data/SingleMuon/UL2018_MiniAODv2_NanoAODv9-v2/2018C/D7033346-A91B-5440-9B0A-8FCB0D90EF4B.root']
        # in_file_names = ['/cms/data/store/user/abrinke1/NanoAOD/2018/data/SingleMuon/UL2018_MiniAODv2_NanoAODv9-v2/2018C/D0D7BD3D-5626-6C4B-B9F4-8DBD9AE223A9.root']
        # in_file_names = ['/cms/data/store/user/abrinke1/NanoAOD/2018/data/SingleMuon/UL2018_MiniAODv2_NanoAODv9-v2/2018C/2B07B4C0-852B-9B4F-83FA-CA6B047542D1.root']
        iFile = 1

    ## TODO: should investigate options - AWB 2022.09.19
    runner = processor.Runner(
        executor=processor.FuturesExecutor(workers=4),
        schema=nanoevents.NanoAODSchema,
        # chunksize=2 ** 19,  ## Originial setting
        # chunksize=2 ** 20,  ## Governs the number of times LeptonJetProcessor "process" is called
        chunksize=2 ** 21,  ## Governs the number of times LeptonJetProcessor "process" is called
    )

    output = runner(
        fileset={LABEL: in_file_names},
        treename="Events",
        processor_instance=LeptonJetProcessor(),
    )

    ## Concatenate histograms from various processes
    out_hists = {}
    for out_h in output:
        h_name = out_h.axes[0].name
        if not h_name in out_hists.keys():
            if VERBOSE: print('  * Creating new histogram %s' % h_name)
            out_hists[h_name] = out_h
        else:
            if VERBOSE: print('  * Adding histogram %s to %s' % (h_name, out_hists[h_name].axes[0].name))
            out_hists[h_name] += out_h

    ## Create output directories and output ROOT file
    if TEST: LABEL += '_TEST'
    #KH out_root = uproot3.recreate('plots/LeptonJetC_%s_%s.root' % (CAT, LABEL))
    out_root = uproot.recreate('plots/LeptonJetC_%s_%s.root' % (CAT, LABEL))
    if not os.path.exists('plots/png/LeptonJetC_%s_%s/' % (CAT, LABEL)):
        os.makedirs('plots/png/LeptonJetC_%s_%s/' % (CAT, LABEL))

    ## Write output pkl file, ROOT file, and png histograms
    print('\nCollected all histograms (?), writing to plots and LeptonJetC_%s_%s.pkl' % (CAT, LABEL))
    with open('plots/LeptonJetC_%s_%s.pkl' % (CAT, LABEL), 'wb') as out_file:
        for h_name in out_hists.keys():
            pickle.dump(out_hists[h_name], out_file)
            print('  * Plotting histogram %s' % h_name)
            ## Export to ROOT file
            # out_root[h_name] = coffea_hist.export1d(out_hists[h_name])
            out_root[h_name] = out_hists[h_name].to_numpy()
            ## Plot to png
            fig = plt.figure()
            out_hists[h_name].plot()
            fig.savefig('plots/png/LeptonJetC_%s_%s/%s.png' % (CAT, LABEL, h_name))
            fig.clear()
            plt.close(fig)
        ## End loop: for h_name in out_hists.keys()
    ## End conditional: with open('plots/LeptonJetC_%s_%s.pkl' % (CAT, LABEL), 'wb') as out_file
    out_root.close()

    print('All done!!!\n')
