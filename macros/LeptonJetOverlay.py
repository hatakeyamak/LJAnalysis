#! /usr/bin/env python

##############################################
###  Overlay histograms from LeptonJet.py  ###
##############################################

import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

import os
import math

#CAT   = 'SingleMu_twoLep_hiM'
#LABEL = '_1000M'
CAT   = 'SingleMu_oneLep_lepMVA0p8_pMET30_eq1j'
LABEL = '_1000M'
DATA  = True


def main():

    print '\nInside LeptonJetOverlay.py\n'

    ## Open input files with histograms
    if DATA: file1_str = 'plots/LeptonJet_%s_SingleMuon_2018C%s.root' % (CAT, LABEL)
    file2_str = 'plots/LeptonJet_%s_QCD_MuEnriched%s.root' % (CAT, LABEL)
    file3_str = 'plots/LeptonJet_%s_WJetsToLNu%s.root' % (CAT, LABEL)
    file4_str = 'plots/LeptonJet_%s_DYJetsToLL_M-10to50%s.root' % (CAT, LABEL)
    file5_str = 'plots/LeptonJet_%s_DYJetsToLL_M-50%s.root' % (CAT, LABEL)
    file6_str = 'plots/LeptonJet_%s_LQsc-SM_ktMin_20%s.root' % (CAT, LABEL)
    # file6_str = 'plots/LeptonJet_%s_LQsc-SM_up_mu_ktMin_20%s.root' % (CAT, LABEL)
    # file7_str = 'plots/LeptonJet_%s_LQsc-SM_up_mu_ktMin_60%s.root' % (CAT, LABEL)

    if DATA:
        print 'Opening file: %s' % file1_str
        file1 = R.TFile(file1_str, 'read')
    print 'Opening file: %s\n' % file2_str
    file2 = R.TFile(file2_str, 'read')
    print 'Opening file: %s\n' % file3_str
    file3 = R.TFile(file3_str, 'read')
    print 'Opening file: %s\n' % file4_str
    file4 = R.TFile(file4_str, 'read')
    print 'Opening file: %s\n' % file5_str
    file5 = R.TFile(file5_str, 'read')
    print 'Opening file: %s\n' % file6_str
    file6 = R.TFile(file6_str, 'read')
    # print 'Opening file: %s\n' % file7_str
    # file7 = R.TFile(file7_str, 'read')

    ## Make directory for output plots
    out_dir = 'plots/LeptonJet_%s%s/png/' % (CAT, LABEL)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## Get list of histograms
    hists = []
    if DATA: file1.cd()
    else:    file2.cd()
    for key in R.gDirectory.GetListOfKeys():
        hists.append(key.GetName())

    ## Make a canvas to draw to
    canv = R.TCanvas('canv')
    canv.cd()

    ## Loop over histograms and overlay
    for hist in hists:

        ## Get histograms
        try:
            if DATA: h1 = file1.Get(hist).Clone()
            h2 = file2.Get(hist).Clone()
            h3 = file3.Get(hist).Clone()
            h4 = file4.Get(hist).Clone()
            h5 = file5.Get(hist).Clone()
            h6 = file6.Get(hist).Clone()
            # h7 = file7.Get(hist).Clone()
        except:
            print '\n*** BIZZARE CASE!!! %s missing from at least one file! ***\n' % hist
            continue

        ## Scale QCD contribution to match data normalization
        # data_sum      = h1.Integral()
        # prompt_MC_sum = h3.Integral() + h4.Integral() + h5.Integral()
        # print 'In histogram %s, scale by %.4f' % (hist, (data_sum - prompt_MC_sum) / h2.Integral())
        h2.Scale(0.5)

        ## Add histograms into a stack
        h3 = h2 + h3
        h4 = h3 + h4
        h5 = h4 + h5

        if h5.Integral() == 0: continue

        ## Scale signal to MC normalization
        print '\nExpect a total of %.2f muon + jet (kT > 20) signal events' % h6.Integral()
        h6.Scale( h5.Integral() / h6.Integral() )
        # h7.Scale( h5.Integral() / h7.Integral() )


        ## Figure out re-binning factor
        ## ## Require an average of > 20 events per filled bin
        nBinsFilled = 0
        for i in range(1, h5.GetNbinsX()+1):
            ## Only count bins as "filled" if they have at least 1% of the average occupancy
            nBinsFilled += 1*(h5.GetBinContent(i) > 0.01*(h5.Integral() / h5.GetNbinsX()))
        rebin = 1
        nIter = 0
        # while ( ((h5.Integral()*rebin / nBinsFilled) < 20) or ((h5.GetNbinsX() / rebin) > 50) ):
        # while (h5.GetNbinsX() / rebin) > 49:
        while (nBinsFilled / rebin) > 101:
        # while (nBinsFilled / rebin) > math.sqrt(h5.Integral()) / 100:
            nIter += 1
            if   (h5.GetNbinsX() / rebin) % 2 == 0:
                rebin *= 2
            elif (h5.GetNbinsX() / rebin) % 3 == 0:
                rebin *= 3
            elif (h5.GetNbinsX() / rebin) % 5 == 0:
                rebin *= 5
            if nIter >= 10: break
        if DATA: h1.Rebin(rebin)
        h2.Rebin(rebin)
        h3.Rebin(rebin)
        h4.Rebin(rebin)
        h5.Rebin(rebin)
        h6.Rebin(rebin)
        # h7.Rebin(rebin)

        # ## Scale data to MC
        # if DATA: h1.Scale( h5.Integral() / h1.Integral() )

        if DATA: h1.SetLineWidth(3)
        if DATA: h1.SetLineColor(R.kBlack)
        h2.SetLineColor(R.kViolet)
        h3.SetLineColor(R.kGreen)
        h4.SetLineColor(R.kAzure-2)
        h5.SetLineColor(R.kAzure+8)

        h6.SetLineColor(R.kRed)
        # h7.SetLineColor(R.kMagenta)
        h6.SetLineWidth(3)
        # h7.SetLineWidth(3)

        h2.SetFillColor(R.kViolet)
        h3.SetFillColor(R.kGreen)
        h4.SetFillColor(R.kAzure-2)
        h5.SetFillColor(R.kAzure+8)

        ## Truncate x-axis range to significant events (> 0.5% of the maximum bin occupancy)
        iMin = 1
        iMax = h5.GetNbinsX()
        for iBin in range(1, iMax):
            if h5.GetBinContent(iBin) > 0.005*h5.GetMaximum() or h6.GetBinContent(iBin) > 0.005*h6.GetMaximum():
                iMin = iBin
                break
        for iBin in reversed(range(iMin+1, iMax+1)):
            if h5.GetBinContent(iBin) > 0.005*h5.GetMaximum() or h6.GetBinContent(iBin) > 0.005*h6.GetMaximum():
                iMax = iBin
                break
        h5.GetXaxis().SetRangeUser( h5.GetBinLowEdge(iMin), h5.GetBinLowEdge(iMax+1) )
        h5.GetYaxis().SetRangeUser( 0, max(h5.GetMaximum(), h6.GetMaximum())*1.2)
        if DATA: h5.GetYaxis().SetRangeUser( 0, max(h1.GetMaximum(), h5.GetMaximum())*1.2)
        h5.Draw('hist')
        h4.Draw('histsame')
        h3.Draw('histsame')
        h2.Draw('histsame')
        if DATA: h1.Draw('histesame')
        h6.Draw('histsame')
        # h7.Draw('histsame')

        canv.SaveAs(out_dir+hist+'.png')


        # ## Compute S^2 / B for Z+jets vs. W+jets
        # ZJets = h5 - h3
        # WJets = h3 - h2
        # ZJets.Scale(1.0 / ZJets.Integral())
        # WJets.Scale(1.0 / WJets.Integral())
        # SigSq2 = 0
        # SigSq20 = 0
        # SigSq200 = 0
        # for i in range(1, h1.GetNbinsX()+1):
        #     # if h1.GetBinContent(i) * h1.GetEntries() > 16:
        #     if WJets.GetBinContent(i) > 0.005:
        #         SigSq2 += pow( ZJets.GetBinContent(i), 2) / WJets.GetBinContent(i)
        #     # if h1.GetBinContent(i) * h1.GetEntries() > 100:
        #     if WJets.GetBinContent(i) > 0.01:
        #         SigSq20 += pow( ZJets.GetBinContent(i), 2) / WJets.GetBinContent(i)
        #     if WJets.GetBinContent(i) > 0.02:
        #         SigSq200 += pow( ZJets.GetBinContent(i), 2) / WJets.GetBinContent(i)
        # if 'MET' in hist:
        #     print '\nFor %s, SigSq2 = %.3f (%.3f, %.3f)' % (hist, SigSq2, SigSq20, SigSq200)
        # del ZJets
        # del WJets


        if DATA: del h1
        del h2
        del h3
        del h4
        del h5
        del h6
        # del h7

    ## Close input files
    if DATA: file1.Close()
    file2.Close()
    file3.Close()
    file4.Close()
    file5.Close()
    file6.Close()
    # file7.Close()

    print '\nAll done!'


## Define 'main' function as primary executable
if __name__ == '__main__':
    main()
