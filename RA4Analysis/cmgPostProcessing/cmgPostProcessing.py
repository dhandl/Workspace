import ROOT
import pickle
import sys, os, copy, random, subprocess, datetime
from array import array
from Workspace.RA4Analysis.cmgObjectSelection import cmgLooseLepIndices, splitIndList, get_cmg_jets_fromStruct, splitListOfObjects, cmgTightMuID, cmgTightEleID
from Workspace.HEPHYPythonTools.xsec import xsec
from Workspace.HEPHYPythonTools.helpers import getObjFromFile, getObjDict, getFileList
from Workspace.HEPHYPythonTools.convertHelpers import compileClass, readVar, printHeader, typeStr, createClassString

from math import *
from Workspace.HEPHYPythonTools.user import username

ROOT.gROOT.ProcessLine(".L ../../HEPHYPythonTools/scripts/root/WPolarizationVariation.C+")
ROOT.gROOT.ProcessLine(".L ../../HEPHYPythonTools/scripts/root/TTbarPolarization.C+")
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()

from Workspace.HEPHYPythonTools.helpers import getChunks
from Workspace.RA4Analysis.cmgTuples_Data25ns_miniAODv2 import *
from Workspace.RA4Analysis.cmgTuples_Spring15_MiniAODv2_25ns import *
from systematics_helper import calc_btag_systematics, calc_LeptonScale_factors_and_systematics
from btagEfficiency import *

bTagEffFile = '/data/dspitzbart/Results2015/MCEff_skim_pkl'

try:
  mcEffDict = pickle.load(file(bTagEffFile))
except IOError:
  print 'Unable to load MC efficiency file!'
  mcEffDict = False

#maxConsideredBTagWeight = 2
#calcSystematics = True
separateBTagWeights = True

defSampleStr = "TTJets_25ns"

subDir = "postProcessing_Syst"

#branches to be kept for data and MC
branchKeepStrings_DATAMC = ["run", "lumi", "evt", "isData", "rho", "nVert",
                     "nJet25", "nBJetLoose25", "nBJetMedium25", "nBJetTight25", "nJet40", "nJet40a", "nBJetLoose40", "nBJetMedium40", "nBJetTight40", 
                     "nLepGood20", "nLepGood15", "nLepGood10", "htJet25", "mhtJet25", "htJet40j", "htJet40", "mhtJet40", "nSoftBJetLoose25", "nSoftBJetMedium25", "nSoftBJetTight25", 
                     "met*","Flag_*","HLT_*",
#                     "nFatJet","FatJet_*", 
                     "nJet", "Jet_*", 
                     "nLepGood", "LepGood_*", 
                     "nLepOther", "LepOther_*", 
                     "nTauGood", "TauGood_*",
                     ] 

#branches to be kept for MC samples only
branchKeepStrings_MC = [ "nTrueInt", "genWeight", "xsec", "puWeight", 
                     "GenSusyMScan1", "GenSusyMScan2", "GenSusyMScan3", "GenSusyMScan4", "GenSusyMGluino", "GenSusyMGravitino", "GenSusyMStop", "GenSusyMSbottom", "GenSusyMStop2", "GenSusyMSbottom2", "GenSusyMSquark", "GenSusyMNeutralino", "GenSusyMNeutralino2", "GenSusyMNeutralino3", "GenSusyMNeutralino4", "GenSusyMChargino", "GenSusyMChargino2", 
                     "ngenLep", "genLep_*", 
                     "nGenPart", "GenPart_*",
                     "ngenPartAll","genPartAll_*" ,
                     "ngenTau", "genTau_*", 
                     "ngenLepFromTau", "genLepFromTau_*"]

#branches to be kept for data only
branchKeepStrings_DATA = []

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--samples", dest="allsamples", default=defSampleStr, type="string", action="store", help="samples:Which samples.")
parser.add_option("--inputTreeName", dest="inputTreeName", default="treeProducerSusySingleLepton", type="string", action="store", help="samples:Which samples.")
parser.add_option("--targetDir", dest="targetDir", default="/data/"+username+"/cmgTuples/"+subDir+'/', type="string", action="store", help="target directory.")
parser.add_option("--skim", dest="skim", default="", type="string", action="store", help="any skim condition?")
parser.add_option("--small", dest="small", default = False, action="store_true", help="Just do a small subset.")
parser.add_option("--overwrite", dest="overwrite", default = False, action="store_true", help="Overwrite?")
parser.add_option("--calcbtagweights", dest="systematics", default = False, action="store_true", help="Calculate b-tag weights for systematics?")
parser.add_option("--btagWeight", dest="btagWeight", default = 2, action="store", help="Max nBJet to calculate the b-tag weight for")
parser.add_option("--hadronicLeg", dest="hadronicLeg", default = False, action="store_true", help="Use only the hadronic leg of the sample?")
parser.add_option("--manScaleFactor", dest="manScaleFactor", default = 1, action="store", help="define a scale factor for the whole sample")

(options, args) = parser.parse_args()
skimCond = "(1)"
ht500lt250 = "Sum$(Jet_pt)>500&&(LepGood_pt[0]+met_pt)>250"
common_skim = "HT500LT250"
if options.skim.startswith('met'):
  skimCond = "met_pt>"+str(float(options.skim[3:]))
if options.skim=='HT400':
  skimCond = "Sum$(Jet_pt)>400"
if options.skim=='HT400ST200':   ##tuples have already ST200 skim
  skimCond = "Sum$(Jet_pt)>400&&(LepGood_pt[0]+met_pt)>200"
if options.skim=='HT500ST250':  
  skimCond = ht500lt250
if options.skim=='LHEHT600':
  skimCond = "lheHTIncoming<600"

skimCond += "&&Sum$(LepGood_pt>25&&abs(LepGood_eta)<2.5)>=0"

##skim conditions for fancy ttJets combination##

####dilep skim##
if options.skim=='HT500ST250diLep':
  skimCond = "((ngenLep+ngenTau)==2)&&lheHTIncoming<=1000&&"+ht500lt250
###semilep skim###
if options.skim=='HT500ST250semiLep':
  skimCond = "((ngenLep+ngenTau)==1)&&lheHTIncoming<=1000&&"+ht500lt250
###Full hadronic###
if options.skim=='HT500ST250LHE_FullHadronic_inc':
  skimCond = "((ngenLep+ngenTau)==0)&&lheHTIncoming<=600&&"+ht500lt250
###Full hadronic for the ht binned###
if options.skim=='HT500ST250LHE_FullHadronic':
  skimCond = "lheHTIncoming>600&&lheHTIncoming<=1000&&((ngenLep+ngenTau)==0)&&"+ht500lt250
###Full inclusive for high HT
if options.skim=='LHEHT1000':
  skimCond = "lheHTIncoming>1000&&"+ht500lt250


if options.hadronicLeg:
  skimCond += "&&(ngenLep+ngenTau)==0"

if options.manScaleFactor!=1:
  target_lumi = target_lumi*float(options.manScaleFactor)
  print
  print "target lumi scaled!"
  print "New lumi:", target_lumi

if options.skim=='inc':
  skimCond = "(1)"

if sys.argv[0].count('ipython'):
  options.small=True

###For PU reweight###
PU_File = ROOT.TFile("/data/easilar/tuples_from_Artur/JECv6recalibrateMET_2100pb/trig_skim/PUhistos/ratio_PU.root")
PU_histo = PU_File.Get("h_ratio")
#####################
###For Lepton SF#####
mu_mediumID_File = ROOT.TFile("/data/easilar/SF2015/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root")
mu_looseID_File = ROOT.TFile("/data/easilar/SF2015/TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta-2.root")
mu_miniIso02_File = ROOT.TFile("/data/easilar/SF2015/TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root")
mu_sip3d_File = ROOT.TFile("/data/easilar/SF2015/TnP_MuonID_NUM_TightIP3D_DENOM_LooseID_VAR_map_pt_eta.root")
ele_kin_File = ROOT.TFile("/data/easilar/SF2015/kinematicBinSFele.root")
#
histos_LS = {
'mu_mediumID_histo':  mu_mediumID_File.Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass"),\
'mu_looseID_histo':   mu_looseID_File.Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass"),\
'mu_miniIso02_histo': mu_miniIso02_File.Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass"),\
'mu_sip3d_histo':     mu_sip3d_File.Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass"),\
'ele_cutbased_histo': ele_kin_File.Get("CutBasedTight"),\
'ele_miniIso01_histo':ele_kin_File.Get("MiniIso0p1_vs_AbsEta"),\
}
#####################

maxConsideredBTagWeight = options.btagWeight
calcSystematics = options.systematics

def getTreeFromChunk(c, skimCond, iSplit, nSplit):
  if not c.has_key('file'):return
  rf = ROOT.TFile.Open(c['file'])
  assert not rf.IsZombie()
  rf.cd()
  tc = rf.Get("tree")
  nTot = tc.GetEntries()
  fromFrac = iSplit/float(nSplit)
  toFrac   = (iSplit+1)/float(nSplit)
  start = int(fromFrac*nTot)
  stop  = int(toFrac*nTot)
  ROOT.gDirectory.cd('PyROOT:/')
  print "Copy tree from source: total number of events found:",nTot,"Split counter: ",iSplit,"<",nSplit,"first Event:",start,"nEvents:",stop-start
  t = tc.CopyTree(skimCond,"",stop-start,start)
  tc.Delete()
  del tc
  rf.Close()
  del rf
  return t
   
def getGenWandLepton(c):
  genPartAll = [getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','motherId','motherIndex'], j) for j in range(int(c.GetLeaf('nGenPart').GetValue()))]
  lepton = filter(lambda l:abs(l['pdgId']) in [11,13,15], genPartAll)
  if len(lepton)==0:
    print "no generated lepton found!"
    p4w=False
    p4lepton=False
    return p4w, p4lepton
  lFromW = filter(lambda w:abs(w['motherId'])==24, lepton)
  if len(lFromW)==0:
    print 'no generated W found!'
    print lepton
    p4w=False
    p4lepton=False
    return p4w, p4lepton
  elif len(lFromW)>0:
    if len(lFromW)>1: print 'this should not have happened'
    if abs(lFromW[0]['motherId'])!=24: print 'this should not have happened'
    genW = getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','motherId','motherIndex'], int(lFromW[0]['motherIndex']))
    lep = ROOT.TLorentzVector()
    lep.SetPtEtaPhiM(lFromW[0]['pt'],lFromW[0]['eta'],lFromW[0]['phi'],lFromW[0]['mass'])
  if abs(genW['pdgId'])!=24: 'this should not have happened'
  W = ROOT.TLorentzVector()
  W.SetPtEtaPhiM(genW['pt'],genW['eta'],genW['phi'],genW['mass'])
  p4lepton = ROOT.LorentzVector(lep.Px(),lep.Py(),lep.Pz(),lep.E())
  p4w = ROOT.LorentzVector(W.Px(),W.Py(),W.Pz(),W.E())
  return p4w, p4lepton

def getGenTopWLepton(c):
  genPartAll = [getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','charge','motherId','motherIndex'], j) for j in range(int(c.GetLeaf('nGenPart').GetValue()))]
  lepton = filter(lambda l:abs(l['pdgId']) in [11,13,15], genPartAll)
  if len(lepton)==0:
    p4t=False
    p4w=False
    p4lepton=False
    return p4t, p4w, p4lepton
  lFromW = filter(lambda w:abs(w['motherId'])==24, lepton)
  if len(lFromW)>0:
    if len(lFromW)==1:
      if abs(lFromW[0]['motherId'])!=24: print '1)this should not have happened'
      genW = getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','charge','motherId','motherIndex'], int(lFromW[0]['motherIndex']))
      if abs(genW['pdgId'])!=24: '2)this should not have happened'
      genTop = getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','charge','motherId','motherIndex'], int(genW['motherIndex']))
      lep = ROOT.TLorentzVector()
      lep.SetPtEtaPhiM(lFromW[0]['pt'],lFromW[0]['eta'],lFromW[0]['phi'],lFromW[0]['mass'])
    elif len(lFromW)==2:
      match = False
      leadLep = getObjDict(c, 'LepGood_', ['pt','eta','phi','mass','pdgId','charge'], 0)
      for l in lFromW:
        if leadLep['charge'] == l['charge']:
          match = True
          genW = getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','charge','motherId','motherIndex'], int(l['motherIndex']))
          genTop = getObjDict(c, 'GenPart_', ['pt','eta','phi','mass','pdgId','charge','motherId','motherIndex'], int(genW['motherIndex']))
          lep = ROOT.TLorentzVector()
          lep.SetPtEtaPhiM(l['pt'],l['eta'],l['phi'],l['mass'])
      if not match:
        print 'No match at all!'
        p4t=False
        p4w=False
        p4lepton=False
        return p4t, p4w, p4lepton
  elif len(lFromW)>2 or len(lFromW)==0:
    print "8) this should not have happened"
    p4t=False
    p4w=False
    p4lepton=False
    return p4t, p4w, p4lepton
  t = ROOT.TLorentzVector()
  W = ROOT.TLorentzVector()
  W.SetPtEtaPhiM(genW['pt'],genW['eta'],genW['phi'],genW['mass'])
  t.SetPtEtaPhiM(genTop['pt'],genTop['eta'],genTop['phi'],genTop['mass'])
  p4lepton = ROOT.LorentzVector(lep.Px(),lep.Py(),lep.Pz(),lep.E())
  p4w = ROOT.LorentzVector(W.Px(),W.Py(),W.Pz(),W.E())
  p4t = ROOT.LorentzVector(t.Px(),t.Py(),t.Pz(),t.E())
  return p4t, p4w, p4lepton

def cleanJetsAndLeptons(jets,leptons,deltaR,arbitration):
    dr2 = deltaR**2
    goodjet = [ True for j in jets ]
    goodlep = [ True for l in leptons ]
    for il, l in enumerate(leptons):
        ibest, d2m = -1, dr2
        for i,j in enumerate(jets):
            d2i = deltaR2(l.eta(),l.phi(), j.eta(),j.phi())
            if d2i < dr2:
                choice = arbitration(j,l)
                if choice == j:
                   # if the two match, and we prefer the jet, then drop the lepton and be done
                   goodlep[il] = False
                   break
                elif choice == (j,l) or choice == (l,j):
                   # asked to keep both, so we don't consider this match
                   continue
            if d2i < d2m:
                ibest, d2m = i, d2i
        # this lepton has been killed by a jet, then we clean the jet that best matches it
        if not goodlep[il]: continue
        if ibest != -1: goodjet[ibest] = False
    return ( [ j for (i ,j) in enumerate(jets)    if goodjet[i ] == True ],
             [ l for (il,l) in enumerate(leptons) if goodlep[il] == True ] )

exec('allSamples=['+options.allsamples+']')
for isample, sample in enumerate(allSamples):
  chunks, sumWeight = getChunks(sample)
  outDir = options.targetDir+"/".join([common_skim, sample['name']])
  if os.path.exists(outDir) and os.listdir(outDir) != [] and not options.overwrite:
    print "Found non-empty directory: %s -> skipping!"%outDir
    continue

  tmpDir = outDir+'/tmp/'
  os.system('mkdir -p ' + outDir) 
  os.system('mkdir -p '+tmpDir)
  os.system('rm -rf '+tmpDir+'/*')
  
  if sample['isData']: 
    lumiScaleFactor=1
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA 
  else:
<<<<<<< HEAD
    if "TTJets" in sample['dbsName']: lumiScaleFactor = xsec[sample['dbsName']]*target_lumi/float(sumWeight)
    else: lumiScaleFactor = target_lumi/float(sumWeight)
    #lumiScaleFactor = target_lumi/float(sumWeight)
=======
    if ("TTJets" in sample['dbsName']): lumiScaleFactor = xsec[sample['dbsName']]*target_lumi/float(sumWeight)
    else: lumiScaleFactor = target_lumi/float(sumWeight)
>>>>>>> origin/74X-master
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC
  
  sampleKey = ''
  if 'TTJets' in sample['dbsName']: sampleKey = 'TTJets'
  elif 'WJets' in sample['dbsName']: sampleKey = 'WJets'
  else: sampleKey = 'none'
  
  readVariables = ['met_pt/F', 'met_phi/F']
  newVariables = ['weight/F','muonDataSet/I','eleDataSet/I']
  aliases = [ "met:met_pt", "metPhi:met_phi"]

  readVectors = [\
<<<<<<< HEAD
    {'prefix':'LepGood', 'nMax':8, 'vars':['pt/F', 'eta/F', 'phi/F', 'pdgId/I', 'relIso03/F','SPRING15_25ns_v1/I' ,'tightId/I', 'miniRelIso/F','mass/F','sip3d/F','mediumMuonId/I', 'mvaIdPhys14/F','mvaIdSpring15/F','lostHits/I', 'convVeto/I', 'charge/I']},
    #{'prefix':'LepGood',  'nMax':8, 'vars':['pt/F', 'eta/F', 'phi/F', 'pdgId/I', 'relIso03/F', 'tightId/I', 'miniRelIso/F','mass/F','sip3d/F','mediumMuonId/I', 'mvaIdPhys14/F','lostHits/I', 'convVeto/I']},
    {'prefix':'Jet',  'nMax':100, 'vars':['pt/F', 'eta/F', 'phi/F', 'id/I','btagCSV/F', 'btagCMVA/F']},
  ]
  if not sample['isData']: 
    newVariables.extend(['weight_XSecTTBar1p1/F','weight_XSecTTBar0p9/F','weight_WPolPlus10/F', 'weight_WPolMinus10/F', 'weight_TTPolPlus5/F', 'weight_TTPolMinus5/F'])
=======
    {'prefix':'LepGood', 'nMax':8, 'vars':['pt/F', 'eta/F', 'phi/F', 'pdgId/I', 'relIso03/F','SPRING15_25ns_v1/I' ,'tightId/I', 'miniRelIso/F','mass/F','sip3d/F','mediumMuonId/I', 'mvaIdPhys14/F','mvaIdSpring15/F','lostHits/I', 'convVeto/I']},
    {'prefix':'Jet',  'nMax':100, 'vars':['pt/F', 'eta/F', 'phi/F', 'id/I','btagCSV/F', 'btagCMVA/F']},
  ]
  if not sample['isData']: 
    newVariables.extend(['puReweight_true/F','puReweight_true_Down/F','puReweight_true_Up/F','weight_diLepTTBar0p5/F','weight_diLepTTBar2p0/F','weight_XSecTTBar1p1/F','weight_XSecTTBar0p9/F','weight_XSecWJets1p1/F','weight_XSecWJets0p9/F'])
    newVariables.extend(['lepton_muSF_looseID/D/1.','lepton_muSF_mediumID/D/1.','lepton_muSF_miniIso02/D/1.','lepton_muSF_sip3d/D/1.','lepton_eleSF_cutbasedID/D/1.','lepton_eleSF_miniIso01/D/1.'])
    newVariables.extend(['lepton_muSF_looseID_err/D/0.','lepton_muSF_mediumID_err/D/0.','lepton_muSF_miniIso02_err/D/0.','lepton_muSF_sip3d_err/D/0.','lepton_eleSF_cutbasedID_err/D/0.','lepton_eleSF_miniIso01_err/D/0.'])
>>>>>>> origin/74X-master
    aliases.extend(['genMet:met_genPt', 'genMetPhi:met_genPhi'])
  newVariables.extend( ['nLooseSoftLeptons/I', 'nLooseHardLeptons/I', 'nTightSoftLeptons/I', 'nTightHardLeptons/I'] )
  newVariables.extend( ['deltaPhi_Wl/F','nBJetMediumCSV30/I','nJet30/I','htJet30j/F','st/F'])
  newVariables.extend( ['leptonPt/F','leptonEt/F','leptonMiniRelIso/F','leptonRelIso03/F' ,\
  'leptonEta/F', 'leptonPhi/F','leptonSPRING15_25ns_v1/I/-2','leptonPdg/I/0', 'leptonInd/I/-1',\
 'leptonMass/F', 'singleMuonic/I', 'singleElectronic/I', 'singleLeptonic/I' ]) #, 'mt2w/F'] )
  if calcSystematics:
    #newVariables.extend( ["weightBTag/F", "weightBTag_SF/F", "weightBTag_SF_b_Up/F", "weightBTag_SF_b_Down/F", "weightBTag_SF_light_Up/F", "weightBTag_SF_light_Down/F"])
    for i in range(maxConsideredBTagWeight+1):
      newVariables.extend( ["weightBTag"+str(i)+"/F", "weightBTag"+str(i)+"_SF/F", "weightBTag"+str(i)+"_SF_b_Up/F", "weightBTag"+str(i)+"_SF_b_Down/F", "weightBTag"+str(i)+"_SF_light_Up/F", "weightBTag"+str(i)+"_SF_light_Down/F"])
      #if i>0:
      newVariables.extend( ["weightBTag"+str(i+1)+"p/F", "weightBTag"+str(i+1)+"p_SF/F", "weightBTag"+str(i+1)+"p_SF_b_Up/F", "weightBTag"+str(i+1)+"p_SF_b_Down/F", "weightBTag"+str(i+1)+"p_SF_light_Up/F", "weightBTag"+str(i+1)+"p_SF_light_Down/F"])
  newVars = [readVar(v, allowRenaming=False, isWritten = True, isRead=False) for v in newVariables]

  
  readVars = [readVar(v, allowRenaming=False, isWritten=False, isRead=True) for v in readVariables]
  for v in readVectors:
    readVars.append(readVar('n'+v['prefix']+'/I', allowRenaming=False, isWritten=False, isRead=True))
    v['vars'] = [readVar(v['prefix']+'_'+vvar, allowRenaming=False, isWritten=False, isRead=True) for vvar in v['vars']]

  printHeader("Compiling class to write")
  writeClassName = "ClassToWrite_"+str(isample)
  writeClassString = createClassString(className=writeClassName, vars= newVars, vectors=[], nameKey = 'stage2Name', typeKey = 'stage2Type')
  s = compileClass(className=writeClassName, classString=writeClassString, tmpDir='/data/'+username+'/tmp/')

  readClassName = "ClassToRead_"+str(isample)
  readClassString = createClassString(className=readClassName, vars=readVars, vectors=readVectors, nameKey = 'stage1Name', typeKey = 'stage1Type', stdVectors=False)
  printHeader("Class to Read")
  r = compileClass(className=readClassName, classString=readClassString, tmpDir='/data/'+username+'/tmp/')

  filesForHadd=[]
  if options.small: chunks=chunks[:1]
  for chunk in chunks:
    sourceFileSize = os.path.getsize(chunk['file'])
    nSplit = 1+int(sourceFileSize/(400*10**6)) #split into 400MB
    if nSplit>1: print "Chunk too large, will split into",nSplit,"of appox 400MB"
    for iSplit in range(nSplit):
      cut = "("+skimCond+")&&("+sample['postProcessingCut']+")" if sample.has_key('postProcessingCut') else skimCond
      t = getTreeFromChunk(chunk, cut, iSplit, nSplit)
      if not t: 
        print "Tree object not found:", t
        continue
      t.SetName("Events")
      nEvents = t.GetEntries()
      for v in newVars:
        v['branch'] = t.Branch(v['stage2Name'], ROOT.AddressOf(s,v['stage2Name']), v['stage2Name']+'/'+v['stage2Type'])
      for v in readVars:
        t.SetBranchAddress(v['stage1Name'], ROOT.AddressOf(r, v['stage1Name']))
      for v in readVectors:
        for var in v['vars']:
          t.SetBranchAddress(var['stage1Name'], ROOT.AddressOf(r, var['stage1Name']))
      for a in aliases:
        t.SetAlias(*(a.split(":")))
      print "File: %s Chunk: %s nEvents: %i (skim: %s) condition: %s lumiScaleFactor: %f"%(chunk['file'],chunk['name'], nEvents, options.skim, skimCond, lumiScaleFactor)
      
      for i in range(nEvents):
        if (i%10000 == 0) and i>0 :
          print i,"/",nEvents  , "name:" , chunk['name']
        s.init()
        r.init()
        t.GetEntry(i)
        genWeight = 1 if sample['isData'] else t.GetLeaf('genWeight').GetValue()
<<<<<<< HEAD
        xsectemp = 1 if sample['isData'] else t.GetLeaf('xsec').GetValue()
        if "TTJets" in sample["name"] : 
          s.weight = lumiScaleFactor*genWeight
        else:
          s.weight = lumiScaleFactor*genWeight*xsectemp
=======
        xsec_branch = 1 if sample['isData'] else t.GetLeaf('xsec').GetValue()
        s.weight = lumiScaleFactor*genWeight
        if sample['isData']:
          if "Muon" in sample['name']:
            s.muonDataSet = True
            s.eleDataSet = False
          if "Electron" in sample['name']:
            s.muonDataSet = False
            s.eleDataSet = True
>>>>>>> origin/74X-master

        nTrueInt = t.GetLeaf('nTrueInt').GetValue()
        
        if not sample['isData']:
          s.muonDataSet = False
          s.eleDataSet = False
          s.weight =xsec_branch*lumiScaleFactor*genWeight
          nTrueInt = t.GetLeaf('nTrueInt').GetValue()
          s.puReweight_true = PU_histo.GetBinContent(PU_histo.FindBin(nTrueInt))
          s.puReweight_true_Down = s.puReweight_true*0.95
          s.puReweight_true_Up = s.puReweight_true*1.05
          ngenLep = t.GetLeaf('ngenLep').GetValue()
          ngenTau = t.GetLeaf('ngenTau').GetValue()
          if ("TTJets" in sample['dbsName']):
            s.weight = lumiScaleFactor*genWeight
            s.weight_XSecTTBar1p1 = s.weight*1.1
            s.weight_XSecTTBar0p9 = s.weight*0.9
            if ngenLep+ngenTau == 2:
              s.weight_diLepTTBar2p0 = s.weight*2.0
              s.weight_diLepTTBar0p5 = s.weight*0.5
            else :
              s.weight_diLepTTBar2p0 = s.weight
              s.weight_diLepTTBar0p5 = s.weight
          else :
            s.weight_XSecTTBar1p1 = s.weight
            s.weight_XSecTTBar0p9 = s.weight
<<<<<<< HEAD
        
        if options.leptonSelection.lower() in ['soft','hard']:
          #get all >=loose lepton indices
          looseLepInd = cmgLooseLepIndices(r) 
          #split into soft and hard leptons
          looseSoftLepInd, looseHardLepInd = splitIndList(r.LepGood_pt, looseLepInd, 25.)
          #select tight soft leptons (no special tight ID for now)
          tightSoftLepInd = looseSoftLepInd #No tight soft selection as of yet 
          #select tight hard leptons (use POG ID)
          ###tightHardLepInd = filter(lambda i:(abs(r.LepGood_pdgId[i])==11 and r.LepGood_relIso03[i]<0.14 and r.LepGood_tightId[i]>=3) \
          ###                               or (abs(r.LepGood_pdgId[i])==13 and r.LepGood_relIso03[i]<0.12 and r.LepGood_tightId[i]), looseHardLepInd)
          tightHardLepInd = filter(lambda i:(abs(r.LepGood_pdgId[i])==11 and cmgTightEleID(r,i)) \
                                         or (abs(r.LepGood_pdgId[i])==13 and cmgTightMuID(r,i)), looseHardLepInd)  


          #print "s lepgood pt: " ,s.LepGood_pt[0]
          s.nLooseSoftLeptons = len(looseSoftLepInd)
          s.nLooseHardLeptons = len(looseHardLepInd)
          s.nTightSoftLeptons = len(tightSoftLepInd)
          s.nTightHardLeptons = len(tightHardLepInd)
          #print "tightHardLepInd:" , tightHardLepInd
          vars = ['pt', 'eta', 'phi', 'miniRelIso','relIso03', 'pdgId', 'SPRING15_25ns_v1']
          allLeptons = [getObjDict(t, 'LepGood_', vars, i) for i in looseLepInd]
          looseSoftLep = [getObjDict(t, 'LepGood_', vars, i) for i in looseSoftLepInd] 
          looseHardLep = [getObjDict(t, 'LepGood_', vars, i) for i in looseHardLepInd]
          tightSoftLep = [getObjDict(t, 'LepGood_', vars, i) for i in tightSoftLepInd]
          tightHardLep =  [getObjDict(t, 'LepGood_', vars, i) for i in tightHardLepInd]
          #print "tightHardLep" , tightHardLep 
          leadingLepInd = None
        if options.leptonSelection=='hard':
          if s.nTightHardLeptons>=1:
            leadingLepInd = tightHardLepInd[0]
            #print "highest pt: " , r.LepGood_pt[0]
            s.leptonPt  = r.LepGood_pt[leadingLepInd]
            s.leptonMiniRelIso = r.LepGood_miniRelIso[leadingLepInd]
            s.leptonRelIso03 = r.LepGood_relIso03[leadingLepInd]
            #print s.leptonMiniRelIso ,s.leptonPt, 'met:', r.met_pt, r.nLepGood, r.LepGood_pt[leadingLepInd],r.LepGood_eta[leadingLepInd], r.LepGood_phi[leadingLepInd] , r.LepGood_pdgId[leadingLepInd], r.LepGood_relIso03[leadingLepInd], r.LepGood_tightId[leadingLepInd], r.LepGood_mass[leadingLepInd]
            s.leptonInd = leadingLepInd 
            s.leptonEta = r.LepGood_eta[leadingLepInd]
            s.leptonPhi = r.LepGood_phi[leadingLepInd]
            s.leptonPdg = r.LepGood_pdgId[leadingLepInd]
            s.leptonMass= r.LepGood_mass[leadingLepInd]
            s.leptonSPRING15_25ns_v1= r.LepGood_SPRING15_25ns_v1[leadingLepInd]
            s.st = r.met_pt + s.leptonPt
          s.singleLeptonic = s.nTightHardLeptons==1
          if s.singleLeptonic:
            s.singleMuonic      =  abs(s.leptonPdg)==13
            s.singleElectronic  =  abs(s.leptonPdg)==11
            if "TTJets" in sample["name"]: #W polarization in TTbar
              p4t, p4w, p4lepton = getGenTopWLepton(t)
              if not p4t and not p4w and not p4lepton:
                s.weight = s.weight
                s.weight_WPolPlus10 = s.weight
                s.weight_WPolMinus10 = s.weight
                s.weight_TTPolPlus5 = s.weight
                s.weight_TTPolMinus5 = s.weight
              else:
                cosTheta = ROOT.ttbarPolarizationAngle(p4t, p4w, p4lepton)
                s.weight = s.weight
                s.weight_WPolPlus10 = s.weight
                s.weight_WPolMinus10 = s.weight
                s.weight_TTPolPlus5 = s.weight * (1. + 0.05*(1.-cosTheta)**2) * 1./(1.+0.05*2./3.) * (1./1.0323239521945559)
                s.weight_TTPolMinus5 = s.weight * (1. - 0.05*(1.-cosTheta)**2) * 1./(1.-0.05*2./3.) * (1.034553190276963956)
            elif "WJets" in sample["name"] and not "TTW" in sample["name"]: #W polarization in W+jets
              p4w, p4lepton = getGenWandLepton(t)
              if not p4w and not p4lepton: 
                s.weight = s.weight
                s.weight_WPolPlus10 = s.weight
                s.weight_WPolMinus10 = s.weight
                s.weight_TTPolPlus5 = s.weight
                s.weight_TTPolMinus5 = s.weight
              else:
                cosTheta = ROOT.WjetPolarizationAngle(p4w, p4lepton)
                s.weight = s.weight
                s.weight_WPolPlus10 = s.weight * (1. + 0.1*(1.-cosTheta)**2) * 1./(1.+0.1*2./3.) * (1./1.04923678332724659) 
                s.weight_WPolMinus10 = s.weight * (1. - 0.1*(1.-cosTheta)**2) * 1./(1.-0.1*2./3.) * (1.05627060952003952)
                s.weight_TTPolPlus5 = s.weight
                s.weight_TTPolMinus5 = s.weight 
            else:
              s.weight = s.weight
              s.weight_WPolPlus10 = s.weight
              s.weight_WPolMinus10 = s.weight
              s.weight_TTPolPlus5 = s.weight
              s.weight_TTPolMinus5 = s.weight
          else:
            s.weight = s.weight
            s.weight_WPolPlus10 = s.weight
            s.weight_WPolMinus10 = s.weight
            s.weight_TTPolPlus5 = s.weight
            s.weight_TTPolMinus5 = s.weight

            s.singleMuonic      = False 
            s.singleElectronic  = False 

        if options.leptonSelection=='soft':
          #Select hardest tight lepton among soft leptons
          if s.nTightSoftLeptons>=1:
            leadingLepInd = tightSoftLepInd[0]
  #          print s.leptonPt, r.LepGood_pt[leadingLepInd],r.LepGood_eta[leadingLepInd], leadingLepInd
            s.leptonPt  = r.LepGood_pt[leadingLepInd]
            s.leptonInd = leadingLepInd 
            s.leptonEta = r.LepGood_eta[leadingLepInd]
            s.leptonPhi = r.LepGood_phi[leadingLepInd]
            s.leptonPdg = r.LepGood_pdgId[leadingLepInd]
            s.leptonMass= r.LepGood_mass[leadingLepInd]
            s.st = r.met_pt + s.leptonPt
          s.singleLeptonic = s.nTightSoftLeptons==1
          if s.singleLeptonic:
            s.singleMuonic      =  abs(s.leptonPdg)==13
            s.singleElectronic  =  abs(s.leptonPdg)==11
          else:
            s.singleMuonic      = False 
            s.singleElectronic  = False 
  #      print "Selected",s.leptonPt
        if options.leptonSelection in ['soft','hard']:
          j_list=['eta','pt','phi','btagCMVA', 'btagCSV', 'id']
          #if not sample['isData']: j_list.extend('partonId')
          jets = filter(lambda j:j['pt']>30 and abs(j['eta'])<2.4 and j['id'], get_cmg_jets_fromStruct(r,j_list))
          #print "jets:" , jets
#          lightJets_, bJetsCMVA = splitListOfObjects('btagCMVA', 0.732, jets) 
          lightJets,  bJetsCSV = splitListOfObjects('btagCSV', 0.890, jets)
          #print "bjetsCMVA:" , bJetsCMVA , "bjetsCSV:" ,  bJetsCSV
          s.htJet30j = sum([x['pt'] for x in jets])
          s.nJet30 = len(jets)
#          s.nBJetMediumCMVA30 = len(bJetsCMVA)
          s.nBJetMediumCSV30 = len(bJetsCSV)
          #print "nbjetsCMVA:" , s.nBJetMediumCMVA30  ,"nbjetsCSV:" ,  s.nBJetMediumCSV30
          #s.mt2w = mt2w.mt2w(met = {'pt':r.met_pt, 'phi':r.met_phi}, l={'pt':s.leptonPt, 'phi':s.leptonPhi, 'eta':s.leptonEta}, ljets=lightJets, bjets=bJetsCSV)
          s.deltaPhi_Wl = acos((s.leptonPt+r.met_pt*cos(s.leptonPhi-r.met_phi))/sqrt(s.leptonPt**2+r.met_pt**2+2*r.met_pt*s.leptonPt*cos(s.leptonPhi-r.met_phi))) 
          #print "deltaPhi:" , s.deltaPhi_Wl
  #          print "Warning -> Why can't I compute mt2w?", s.mt2w, len(jets), len(bJets), len(allTightLeptons),lightJets,bJets, {'pt':s.type1phiMet, 'phi':s.type1phiMetphi}, {'pt':s.leptonPt, 'phi':s.leptonPhi, 'eta':s.leptonEta}

        if calcSystematics:
#          separateBTagWeights = False
          zeroTagWeight = 1.
          mceff = getMCEfficiencyForBTagSF(t, mcEffDict[sampleKey], sms='')
          #print
          #print mceff["mceffs"]
          mceffW                = getTagWeightDict(mceff["mceffs"], maxConsideredBTagWeight)
          mceffW_SF             = getTagWeightDict(mceff["mceffs_SF"], maxConsideredBTagWeight)
          mceffW_SF_b_Up        = getTagWeightDict(mceff["mceffs_SF_b_Up"], maxConsideredBTagWeight)
          mceffW_SF_b_Down      = getTagWeightDict(mceff["mceffs_SF_b_Down"], maxConsideredBTagWeight)
          mceffW_SF_light_Up    = getTagWeightDict(mceff["mceffs_SF_light_Up"], maxConsideredBTagWeight)
          mceffW_SF_light_Down  = getTagWeightDict(mceff["mceffs_SF_light_Down"], maxConsideredBTagWeight)
          if not separateBTagWeights:
            lweight = str(s.weight)
          else: lweight = "(1.)"
          #if not separateBTagWeights:
          for i in range(1, maxConsideredBTagWeight+2):
            exec("s.weightBTag"+str(i)+"p="+lweight)
            exec("s.weightBTag"+str(i)+"p_SF="+lweight)
            exec("s.weightBTag"+str(i)+"p_SF_b_Up="+lweight)
            exec("s.weightBTag"+str(i)+"p_SF_b_Down="+lweight)
            exec("s.weightBTag"+str(i)+"p_SF_light_Up="+lweight)
            exec("s.weightBTag"+str(i)+"p_SF_light_Down="+lweight)
          for i in range(maxConsideredBTagWeight+1):
            exec("s.weightBTag"+str(i)+"="              +str(mceffW[i])+'*'+lweight)
            exec("s.weightBTag"+str(i)+"_SF="           +str(mceffW_SF[i])+'*'+lweight)
            exec("s.weightBTag"+str(i)+"_SF_b_Up="      +str(mceffW_SF_b_Up[i])+'*'+lweight)
            exec("s.weightBTag"+str(i)+"_SF_b_Down="    +str(mceffW_SF_b_Down[i])+'*'+lweight)
            exec("s.weightBTag"+str(i)+"_SF_light_Up="  +str(mceffW_SF_light_Up[i])+'*'+lweight)
            exec("s.weightBTag"+str(i)+"_SF_light_Down="+str(mceffW_SF_light_Down[i])+'*'+lweight)
            for j in range(i+1, maxConsideredBTagWeight+2):
              exec("s.weightBTag"+str(j)+"p               -="+str(mceffW[i])+'*'+lweight) #prob. for >=j b-tagged jets
              exec("s.weightBTag"+str(j)+"p_SF            -="+str(mceffW_SF[i])+'*'+lweight)
              exec("s.weightBTag"+str(j)+"p_SF_b_Up       -="+str(mceffW_SF_b_Up[i])+'*'+lweight)
              exec("s.weightBTag"+str(j)+"p_SF_b_Down     -="+str(mceffW_SF_b_Down[i])+'*'+lweight)
              exec("s.weightBTag"+str(j)+"p_SF_light_Up   -="+str(mceffW_SF_light_Up[i])+'*'+lweight)
              exec("s.weightBTag"+str(j)+"p_SF_light_Down -="+str(mceffW_SF_light_Down[i])+'*'+lweight)
          for i in range (int(r.nJet)+1, maxConsideredBTagWeight+1):
            exec("s.weightBTag"+str(i)+"               = 0.")
            exec("s.weightBTag"+str(i)+"_SF            = 0.")
            exec("s.weightBTag"+str(i)+"_SF_b_Up       = 0.")
            exec("s.weightBTag"+str(i)+"_SF_b_Down     = 0.")
            exec("s.weightBTag"+str(i)+"_SF_light_Up   = 0.")
            exec("s.weightBTag"+str(i)+"_SF_light_Down = 0.")
            exec("s.weightBTag"+str(i)+"p              = 0.")
            exec("s.weightBTag"+str(i)+"p_SF           = 0.")
            exec("s.weightBTag"+str(i)+"p_SF_b_Up      = 0.")
            exec("s.weightBTag"+str(i)+"p_SF_b_Down    = 0.")
            exec("s.weightBTag"+str(i)+"p_SF_light_Up  = 0.")
            exec("s.weightBTag"+str(i)+"p_SF_light_Down= 0.")
=======
            s.weight_diLepTTBar2p0 = s.weight
            s.weight_diLepTTBar0p5 = s.weight
          if "WJets" in sample['dbsName']:
            s.weight_XSecWJets1p1 = s.weight*1.1
            s.weight_XSecWJets0p9 = s.weight*0.9
          else :
            s.weight_XSecWJets1p1 = s.weight
            s.weight_XSecWJets0p9 = s.weight       
 
        #get all >=loose lepton indices
        looseLepInd = cmgLooseLepIndices(r) 
        #split into soft and hard leptons
        looseSoftLepInd, looseHardLepInd = splitIndList(r.LepGood_pt, looseLepInd, 25.)
        #select tight soft leptons (no special tight ID for now)
        tightSoftLepInd = looseSoftLepInd #No tight soft selection as of yet 
        #select tight hard leptons (use POG ID)
        tightHardLepInd = filter(lambda i:(abs(r.LepGood_pdgId[i])==11 and cmgTightEleID(r,i)) \
                                       or (abs(r.LepGood_pdgId[i])==13 and cmgTightMuID(r,i)), looseHardLepInd)  


        s.nLooseSoftLeptons = len(looseSoftLepInd)
        s.nLooseHardLeptons = len(looseHardLepInd)
        s.nTightSoftLeptons = len(tightSoftLepInd)
        s.nTightHardLeptons = len(tightHardLepInd)
        vars = ['pt', 'eta', 'phi', 'miniRelIso','relIso03', 'pdgId', 'SPRING15_25ns_v1']
        allLeptons = [getObjDict(t, 'LepGood_', vars, i) for i in looseLepInd]
        looseSoftLep = [getObjDict(t, 'LepGood_', vars, i) for i in looseSoftLepInd] 
        looseHardLep = [getObjDict(t, 'LepGood_', vars, i) for i in looseHardLepInd]
        tightSoftLep = [getObjDict(t, 'LepGood_', vars, i) for i in tightSoftLepInd]
        tightHardLep =  [getObjDict(t, 'LepGood_', vars, i) for i in tightHardLepInd]
        leadingLepInd = None
        if s.nTightHardLeptons>=1:
          leadingLepInd = tightHardLepInd[0]
          s.leptonPt  = r.LepGood_pt[leadingLepInd]
          s.leptonMiniRelIso = r.LepGood_miniRelIso[leadingLepInd]
          s.leptonRelIso03 = r.LepGood_relIso03[leadingLepInd]
          s.leptonInd = leadingLepInd 
          s.leptonEta = r.LepGood_eta[leadingLepInd]
          s.leptonPhi = r.LepGood_phi[leadingLepInd]
          s.leptonPdg = r.LepGood_pdgId[leadingLepInd]
          s.leptonMass= r.LepGood_mass[leadingLepInd]
          s.leptonSPRING15_25ns_v1= r.LepGood_SPRING15_25ns_v1[leadingLepInd]
          s.st = r.met_pt + s.leptonPt
        s.singleLeptonic = s.nTightHardLeptons==1
        if s.singleLeptonic:
          lep_vec = ROOT.TLorentzVector()
          lep_vec.SetPtEtaPhiM(s.leptonPt,s.leptonEta,s.leptonPhi,s.leptonMass)
          s.leptonEt = lep_vec.Et()
          s.singleMuonic      =  abs(s.leptonPdg)==13
          s.singleElectronic  =  abs(s.leptonPdg)==11
        else:
          s.singleMuonic      = False 
          s.singleElectronic  = False 

        j_list=['eta','pt','phi','btagCMVA', 'btagCSV', 'id']
        jets = filter(lambda j:j['pt']>30 and abs(j['eta'])<2.4 and j['id'], get_cmg_jets_fromStruct(r,j_list))
        lightJets,  bJetsCSV = splitListOfObjects('btagCSV', 0.890, jets)
        s.htJet30j = sum([x['pt'] for x in jets])
        s.nJet30 = len(jets)
        s.nBJetMediumCSV30 = len(bJetsCSV)
        #s.mt2w = mt2w.mt2w(met = {'pt':r.met_pt, 'phi':r.met_phi}, l={'pt':s.leptonPt, 'phi':s.leptonPhi, 'eta':s.leptonEta}, ljets=lightJets, bjets=bJetsCSV)
        s.deltaPhi_Wl = acos((s.leptonPt+r.met_pt*cos(s.leptonPhi-r.met_phi))/sqrt(s.leptonPt**2+r.met_pt**2+2*r.met_pt*s.leptonPt*cos(s.leptonPhi-r.met_phi))) 


        calc_LeptonScale_factors_and_systematics(s,histos_LS)
        if calcSystematics: 
          calc_btag_systematics(t,s,r,mcEffDict,sampleKey,maxConsideredBTagWeight,separateBTagWeights)
>>>>>>> origin/74X-master

        for v in newVars:
          v['branch'].Fill()
      newFileName = sample['name']+'_'+chunk['name']+'_'+str(iSplit)+'.root'
      filesForHadd.append(newFileName)
      if not options.small:
      #if options.small:
        f = ROOT.TFile(tmpDir+'/'+newFileName, 'recreate')
        t.SetBranchStatus("*",0)
        for b in branchKeepStrings + [v['stage2Name'] for v in newVars] +  [v.split(':')[1] for v in aliases]:
          t.SetBranchStatus(b, 1)
        t2 = t.CloneTree()
        t2.Write()
        f.Close()
        print "Written",tmpDir+'/'+newFileName
        del f
        del t2
        t.Delete()
        del t
      for v in newVars:
        del v['branch']

  print "Event loop end"
  if not options.small: 
    size=0
    counter=0
    files=[]
    for f in filesForHadd:
      size+=os.path.getsize(tmpDir+'/'+f)
      files.append(f)
      if size>(0.5*(10**9)) or f==filesForHadd[-1] or len(files)>200:
        ofile = outDir+'/'+sample['name']+'_'+options.skim+'_'+str(counter)+'.root'
        print "Running hadd on", tmpDir, files
        os.system('cd '+tmpDir+';hadd -f '+ofile+' '+' '.join(files))
        print "Written", ofile
        size=0
        counter+=1
        files=[]
    os.system("rm -rf "+tmpDir)

