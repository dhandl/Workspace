import ROOT
from array import array
from math import *
import os, copy, sys
from Workspace.RA4Analysis.helpers import nameAndCut, nJetBinName, nBTagBinName

ROOT.gROOT.LoadMacro("../../HEPHYPythonTools/scripts/root/tdrstyle.C")
ROOT.setTDRStyle()
ROOT.gStyle.SetOptFit(0000)

ROOT.TH1F().SetDefaultSumw2()
def getRCS(c, cut, dPhiCut):
  h = getPlotFromChain(c, dPhiStr, [0,dPhiCut,pi], cutString=cut, binningIsExplicit=True)
  if h.GetBinContent(1)>0 and h.GetBinContent(2)>0:
    rcs = h.GetBinContent(2)/h.GetBinContent(1)
    rCSE_sim = rcs*sqrt(h.GetBinError(2)**2/h.GetBinContent(2)**2 + h.GetBinError(1)**2/h.GetBinContent(1)**2)
    rCSE_pred = rcs*sqrt(1./h.GetBinContent(2)**2 + 1./h.GetBinContent(1)**2)
    del h
    return {'rCS':rcs, 'rCSE_pred':rCSE_pred, 'rCSE_sim':rCSE_sim}
  del h

from Workspace.HEPHYPythonTools.helpers import getObjFromFile, getChain, getChunks
import Workspace.HEPHYPythonTools.xsec as xsec
#from Workspace.RA4Analysis.simplePlotsCommon import *
from Workspace.RA4Analysis.helpers import *
from Workspace.HEPHYPythonTools.user import username
from draw_helpers import *
from Workspace.RA4Analysis.cmgTuplesPostProcessed_v6_Phys14V2_HT400ST150_withDF import *

ROOT_colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kCyan]

lepSel = 'hard'
#dPhiStr = 'acos((leptonPt+met*cos(leptonPhi-metPhi))/sqrt(leptonPt**2+met**2+2*met*leptonPt*cos(leptonPhi-metPhi)))'
cWJets  = getChain(WJetsHTToLNu[lepSel],histname='')
cTTJets = getChain(ttJets[lepSel],histname='')
cBkg    = getChain([WJetsHTToLNu[lepSel], ttJets[lepSel], QCD[lepSel], DY[lepSel], singleTop[lepSel], TTVH[lepSel]],histname='')

prefix = 'hardSingleLeptonic_ttJetsFit_MET200'
presel='singleLeptonic&&nLooseHardLeptons==1&&nTightHardLeptons==1&&nLooseSoftPt10Leptons==0&&met_pt>200'
#presel = exactlyOneTightLepton(lepton="both")+'&&'+looseLeptonVeto(lepton='both')
uDir = username[0]+'/'+username
subDir = 'pngCMG2/rCS'

#totalYield=0
#for s in samples:
#  cs=getChunks(s, treeName="treeProducerSusySingleLepton")
#  s.update(cs[0][0])
#  s['chain']=ROOT.TChain('tree')
#  s['chain'].Add(s['file'])
#  s["weight"]=4000*xsec.xsec[s['dbsName']]/float(cs[1])
#  nEntry=s['chain'].GetEntries()
#  totalYield+=4000*xsec.xsec[s['dbsName']]/float(cs[1])*nEntry
#  print s["name"],"xsec",xsec.xsec[s['dbsName']],"NSim",cs[1],"nEntry",nEntry, "yield:",4000*xsec.xsec[s['dbsName']]/float(cs[1])*nEntry
#print "totalYield", totalYield

htbins = [(500,750), (750,1000), (1000,1250), (1250,-1)]
#prefix=None
njCorr=0.

res = {}
for lepton in ["both"]:# ["muon", "electron"]: #["both"]:
  if lepton=="muon":
     leptonMinPt=25; minID=1; minRelIso=0.12; leptonMaxEta = 2.4
  if lepton=="electron":
     leptonMinPt=25; minID=3; minRelIso=0.14; leptonMaxEta = 2.4
  if lepton=="both":
     leptonMinPt=(25,25); minID=(3,1); minRelIso=(0.14,0.12); leptonMaxEta = (2.4, 2.4)
  for htb in htbins:
    res[htb] = {}
    for stb in [(250,350), (350,450), (450,-1)]:
      res[htb][stb] = {}
  
      fname = nameAndCut(stb, htb, njetb=None, btb=None)[0]
  #    fname="_".join(([prefix] if prefix else [])+[lepton,name])
      plots={}
      njet_bins = [(4,4),(5,5),(6,-1)]
      nbtag_bins = [0,1,2,3,4]
      for njb in njet_bins: 
        name, cut = nameAndCut(stb, htb, njb, btb=None, presel=presel)
#        cut = "&&".join([
#            exactlyOneTightLepton(lepton=lepton, minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso),\
#            looseLeptonVeto(lepton=lepton, minPt=10), \
#            nJetCut(njb=njb, minPt=30, maxEta=2.4), \
#            htCut  (htb=htb, minPt=30, maxEta=2.4, njCorr=njCorr), \
#            stCut(lepton=lepton, stb=stb, minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso), \
#            ])
        rCS_vs_nbtag = ROOT.TProfile('profile_rCS_nbtag','', len(nbtag_bins)-1, array('d',nbtag_bins), 0, 1)
        print "htb", htb, "stb", stb, "njb", njb, "cut:", cut
        #rCS_vs_nbtag.Reset()
        cTTJets.Draw('(Sum$((deltaPhi_Wl)*('+presel+')>1==1)):Sum$(nBJetMediumCMVA30)>>profile_rCS_nbtag',cut,'goff')
#        ttJets_fromEOS['chain'].Draw(dPhiCut(minDPhi=1, lepton=lepton,minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso)+":"+nBTagStr(minPt=30, maxEta=2.4, minCMVATag=0.732)+'>>profile_rCS_nbtag',cut,'goff') 
#        cTTJets.Draw(dPhiCut(minDPhi=1, lepton=lepton,minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso)+":"+nBTagStr(minPt=30, maxEta=2.4, minCMVATag=0.732)+'>>profile_rCS_nbtag',cut,'goff') 
        plots[njb]=ROOT.gDirectory.Get('profile_rCS_nbtag').Clone()
      opt="eh1"
      c1 = ROOT.TCanvas()
      l=ROOT.TLegend(0.6,1.0-0.07*len(njet_bins),1.0,1.0)
      l.SetFillColor(ROOT.kWhite)
      l.SetShadowColor(ROOT.kWhite)
      l.SetBorderSize(1)
      for i_njb,njb in enumerate(njet_bins):
        ROOT.gStyle.SetOptStat(0)
        plots[njb].GetYaxis().SetRangeUser(0,0.2)
        plots[njb].GetYaxis().SetTitle("R_{CS}")
        plots[njb].GetXaxis().SetTitle("b-tag multiplicity")
        plots[njb].SetLineColor(ROOT_colors[i_njb])
        plots[njb].SetMarkerStyle(0)
        plots[njb].SetMarkerSize(0)
        plots[njb].SetLineStyle(1)
        plots[njb].SetLineWidth(2)
        plots[njb].SetMarkerColor(ROOT_colors[i_njb])
        plots[njb].Draw(opt)
        l.AddEntry(plots[njb], nJetBinName(njb))
        opt="eh1same"
      l.Draw()
#      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/btag_'+fname+'.png')
#      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/btag_'+fname+'.pdf')
#      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/btag_'+fname+'.root')
  
      plots={}
      njet_bins = [4,5,6,7,8,9]
      nbtag_bins = [(0,0), (1,1), (2,2)]
  
      for btb in nbtag_bins: 
        name, cut = nameAndCut(stb,htb,njetb=None,btb=btb, presel=presel, btagVar = 'nBJetMediumCMVA30')
#        cut= "&&".join([
#            exactlyOneTightLepton(lepton=lepton, minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso),\
#            looseLeptonVeto(lepton=lepton, minPt=10), \
#            nBTagCut(btb, minPt=30, maxEta=2.4, minCMVATag=0.732),
#            htCut  (htb=htb, minPt=30, maxEta=2.4, njCorr=njCorr), \
#            stCut(lepton=lepton, stb=stb, minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso), \
#            ])
  
        rCS_vs_njet = ROOT.TProfile('profile_rCS_njet_bTag'+str(btb[0]),'', len(njet_bins)-1, array('d',njet_bins), 0, 1)
        print "htb", htb, "stb", stb, "btb", btb, "cut:", cut
        #rCS_vs_nbtag.Reset()
        cTTJets.Draw('(Sum$((deltaPhi_Wl)*('+presel+')>1==1)):nJet30>>profile_rCS_njet_bTag'+str(btb[0]),cut,'goff')
#        ttJets_fromEOS['chain'].Draw(dPhiCut(minDPhi=1, lepton=lepton,minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso)+":"+nJetStr(minPt=30, maxEta=2.4)+'>>profile_rCS_njet',cut,'goff') 
#        cWJets.Draw(dPhiCut(minDPhi=1, lepton=lepton,minPt=leptonMinPt, maxEta=leptonMaxEta, minID=minID, minRelIso=minRelIso)+":"+nJetStr(minPt=30, maxEta=2.4)+'>>profile_rCS_njet_bTag'+str(btb[0]),cut,'goff') 
        plots[btb]=ROOT.gDirectory.Get('profile_rCS_njet_bTag'+str(btb[0])).Clone()
  
      opt="eh1"
      c1 = ROOT.TCanvas()
      l=ROOT.TLegend(0.6,1.0-0.07*len(nbtag_bins),0.95,0.95)
      l.SetFillColor(ROOT.kWhite)
      l.SetShadowColor(ROOT.kWhite)
      l.SetBorderSize(1)
      for i_btb,btb in enumerate(nbtag_bins):
        rd ={} 
        plots[btb].Fit('pol0','','same')
        FitFunc     = plots[btb].GetFunction('pol0')
        FitPar      = FitFunc.GetParameter(0)
        FitParError = FitFunc.GetParError(0)
        FitFunc.SetLineColor(ROOT_colors[i_btb])
        FitFunc.SetLineStyle(2)
        FitFunc.SetLineWidth(2)
        rd['FitFunction']  = FitFunc
        rd['FitParameter'] = FitPar
        rd['FitParError']  = FitParError
        res[htb][stb][btb] = rd
        ROOT.gStyle.SetOptStat(0)
        plots[btb].GetYaxis().SetRangeUser(0,0.2)
        plots[btb].GetYaxis().SetTitle("R_{CS}")
        plots[btb].GetXaxis().SetTitle("jet multiplicity")
        plots[btb].SetLineColor(ROOT_colors[i_btb])
        plots[btb].SetMarkerStyle(0)
        plots[btb].SetMarkerSize(0)
        plots[btb].SetLineStyle(1)
        plots[btb].SetLineWidth(2)
        plots[btb].SetMarkerColor(ROOT_colors[i_btb])
        plots[btb].Draw(opt)
        l.AddEntry(plots[btb], nBTagBinName(btb))
        opt="eh1same"
  
      l.Draw()
      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/'+prefix+'njet_'+fname+'.png')
      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/'+prefix+'njet_'+fname+'.pdf')
      c1.Print('/afs/hephy.at/user/'+uDir+'/www/'+subDir+'/'+prefix+'njet_'+fname+'.root')

