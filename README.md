# LJAnalysis
Event-loop version
```
/cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_10_6_4 # one time operation
cd CMSSW_10_6_4/src/
cmsenv
git cms-init # one time operation
git clone git@github.com:hatakeyamak/LJAnalysis.git # one time operation
cd LJAnalysis
python2 macros/LeptonJet.py WJetsToLNu SingleMu_oneLep_ge1j 200000 20 1
```

Coffea version
```
# Set up coffea via conda/mamba
git clone git@github.com:hatakeyamak/LJAnalysis.git # one time operation
cd LJAnalysis
python macros/LeptonJetC.py 
```
