

###################################################################################
###  Functions to select events according to common criteria (nPV, JSON, etc.)  ###
###################################################################################

import json

import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn


JSON_FILE = 'data/JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'  ## 2018 Golden JSON

VERBOSE = False

class GetSelection:

    def __init__(self):

        print '\nLoading 2018 Golden JSON from %s' % JSON_FILE
        with open(JSON_FILE) as json_file:
            self.JSON = json.load(json_file)

    ## End def __init__(self)


    ## Determine whether event is in the Golden JSON
    def PassJSON(self, run, LS):

        Good_LS = False

        if not str(run) in self.JSON.keys():
            return Good_LS

        for iLS in range( len( self.JSON[str(run)] ) ):
            if int(LS) >= self.JSON[str(run)][iLS][0] and int(LS) <= self.JSON[str(run)][iLS][-1]:
                Good_LS = True
                break

        return Good_LS
    ## End function: PassJSON(self, run, LS)


## End class GetSelection
