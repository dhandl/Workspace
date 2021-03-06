import ROOT
from ROOT import RooFit as rf
import pickle 
import copy, os, sys
ROOT.gROOT.LoadMacro("../../HEPHYPythonTools/scripts/root/tdrstyle.C")
ROOT.TH1F().SetDefaultSumw2()
ROOT.setTDRStyle()
ROOT.gStyle.SetMarkerStyle(1)
ROOT.gStyle.SetOptTitle(0)

from Workspace.HEPHYPythonTools.helpers import *
from Workspace.HEPHYPythonTools.xsec import *
from Workspace.RA4Analysis.helpers import *
from Workspace.RA4Analysis.signalRegions import *
from Workspace.RA4Analysis.cmgTuples_Spring15_25ns_postProcessed_antiSel import *
#from Workspace.RA4Analysis.cmgTuples_Data_25ns_postProcessed_antiSel import *
from draw_helpers import *
from math import *
from Workspace.HEPHYPythonTools.user import username
from LpTemplateFit import LpTemplateFit

preprefix = 'QCDestimation/final2p1fb/MC'
wwwDir = '/afs/hephy.at/user/'+username[0]+'/'+username+'/www/RunII/Spring15_25ns/'+preprefix+'/'
picklePath = '/data/'+username+'/results2015/QCDEstimation/'
prefix = 'Lp_singleElectronic_'
picklePresel = '20151216_QCDestimation_MC2p1fb_pkl'
pickleFit    = '20151216_fitResult_MC2p1fb_pkl'

if not os.path.exists(wwwDir):
  os.makedirs(wwwDir)

inclusiveTemplate = {(3, 4): {(250,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}}} #use inclusive LT,HT region to get the shape for the fit template

fitCR =  {(3, 4): {(250,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                   (250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}, #QCD CR exclusive in LT and inclusive in HT, where the fits are performed
                   (350,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}, 
                   (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                   (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}}}

signalRegion = {(3, 4): {(250, 350): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}}, #3-4jets W+jets control region
                                      (500, 750):  {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                         (350, -1):  {(500, -1):   {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}}},
                         (350, 450): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0},
                                                    (0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (500, 750):  {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0},
                                                    (0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (500, 750):  {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}}}},
                (4, 5): {(250, 350): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}}, #4-5jets TTbar control region
                                      (500, 750):  {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                         (350, -1):  {(500, -1):   {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}}},
                         (350, 450): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0},
                                                    (0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (500, 750):  {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0},
                                                    (0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (500, 750):  {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'sys':0.025,  'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}}}},
                (5, 5): {(250, 350): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}}},  #signal regions
                         (350, 450): {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'sys':0.025,  'deltaPhi': 1.0},
                                                    (0.75):   {'sys':0.025,  'deltaPhi': 0.75}}}},
                (6, 7): {(250, 350): {(500, 750):  {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                         (350, 450): {(500, 750):  {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.05,   'deltaPhi': 1.0}}},
                          (450, -1): {(500, 750):  {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'sys':0.05,   'deltaPhi': 0.75}}}},
                (8, -1): {(250, 350):{(500, 750):  {(1.0):    {'sys':0.1,   'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'sys':0.1,   'deltaPhi': 1.0}}},
                          (350, -1): {(500, -1):   {(0.75):   {'sys':0.1,   'deltaPhi': 0.75}}},
                          (350, 450):{(500, -1):   {(0.75):   {'sys':0.1,   'deltaPhi': 0.75}}},
                          (450, -1): {(500, -1):   {(0.75):   {'sys':0.1,   'deltaPhi': 0.75}}}}
}

#inclusiveTemplate = {(3, 3): {(250,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}}} #use inclusive LT,HT region to get the shape for the fit template
#
#fitCR =  {(3, 3): {(250,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                   (250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}, #QCD CR exclusive in LT and inclusive in HT, where the fits are performed
#                   #(350,  -1): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}, 
#                   (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                   (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}}}}
#
#signalRegion = {(3, 3): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}, #3-4jets W+jets control region
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 1000): {(0.75):    {'deltaPhi': 0.75}},
#                                      (1000, -1):  {(0.75):    {'deltaPhi': 0.75}}}},
#                (4, 5): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}, #4-5jets TTbar control region
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 1000): {(0.75):    {'deltaPhi': 0.75}},
#                                      (1000, -1):  {(0.75):    {'deltaPhi': 0.75}}}},
#                (4, 4): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
#                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
#                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0}},
#                                      (500, 1000): {(0.75):    {'deltaPhi': 0.75}},
#                                      (1000, -1):  {(0.75):    {'deltaPhi': 0.75}}}}
#}

btreg = [(0,0), (1,1), (2,-1)] #1b and 2b estimates are needed for the btag fit

def makeWeight(lumi=3., sampleLumi=3.,debug=False):
  if debug:
    print 'No lumi-reweighting done!!'
    return 'weight', 'weight*weight'
  else:
    weight_str = '(((weight)/'+str(sampleLumi)+')*'+str(lumi)+')'
    weight_err_str = '('+weight_str+'*'+weight_str+')'
    return weight_str, weight_err_str
lumi = 2.11
sampleLumi = 2.11
debugReweighting = True
weight_str, weight_err_str = makeWeight(lumi, sampleLumi=sampleLumi, debug=debugReweighting)

def getRCS(c, cut, dPhiCut, useWeight = False, weight = 'weight'):
#  dPhiStr = 'acos((LepGood_pt+met_pt*cos(LepGood_phi-met_phi))/sqrt(LepGood_pt**2+met_pt**2+2*met_pt*LepGood_pt*cos(LepGood_phi-met_phi)))'
  dPhiStr = 'deltaPhi_Wl'
  if useWeight:
    h = getPlotFromChain(c, dPhiStr, [0,dPhiCut,pi], cutString=cut, binningIsExplicit=True, weight =  weight)
  else:
    h = getPlotFromChain(c, dPhiStr, [0,dPhiCut,pi], cutString=cut, binningIsExplicit=True, weight='(1)')
  h.Sumw2()
  if h.GetBinContent(1)>0:
    rcs = h.GetBinContent(2)/h.GetBinContent(1)
    if h.GetBinContent(2)>0:
      rCSE_sim = rcs*sqrt(h.GetBinError(2)**2/h.GetBinContent(2)**2 + h.GetBinError(1)**2/h.GetBinContent(1)**2)
      rCSE_pred = rcs*sqrt(1./h.GetBinContent(2) + 1./h.GetBinContent(1))
      del h
      return {'rCS':rcs, 'rCSE_pred':rCSE_pred, 'rCSE_sim':rCSE_sim}
    else:
      del h
      return {'rCS':rcs, 'rCSE_pred':float('nan'), 'rCSE_sim':float('nan')}
  else:
    del h
    return {'rCS':float('nan'), 'rCSE_pred':float('nan'), 'rCSE_sim':float('nan')}

def getPseudoRCS(small,smallE,large,largeE): 
  if small>0:
    rcs = large/small
    if large>0:
      rCSE_sim = rcs*sqrt(smallE**2/small**2 + largeE**2/large**2)
      rCSE_pred = rcs*sqrt(1./small + 1./large)
      return {'rCS':rcs, 'rCSE_pred':rCSE_pred, 'rCSE_sim':rCSE_sim}
    else:
      return {'rCS':rcs, 'rCSE_pred':float('nan'), 'rCSE_sim':float('nan')}
  else:
    return {'rCS':float('nan'), 'rCSE_pred':float('nan'), 'rCSE_sim':float('nan')}

#trigger and filters for real Data
trigger = "&&(HLT_EleHT350||HLT_MuHT350)"
filters = "&&Flag_goodVertices && Flag_HBHENoiseFilter_fix && Flag_CSCTightHaloFilter && Flag_eeBadScFilter && Flag_HBHENoiseIsoFilter"
#filters = "&&Flag_CSCTightHaloFilter&&Flag_HBHENoiseFilter_fix&&Flag_HBHENoiseFilter&&Flag_goodVertices&&Flag_eeBadScFilter&&Flag_EcalDeadCellTriggerPrimitiveFilter"
#filters = "&&Flag_CSCTightHaloFilter&&Flag_HBHENoiseFilter_fix&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_eeBadScFilter"

presel = 'nLep==1&&nVeto==0&&leptonPt>25&&nEl==1&&Jet2_pt>80'
antiSelStr = presel+'&&Selected==(-1)'
SelStr = presel+'&&Selected==1'

cQCD  = getChain(QCDHT_25ns,histname='')
cEWK  = getChain([WJetsHTToLNu_25ns, TTJets_combined_2, singleTop_25ns, DY_25ns, TTV_25ns],histname='')
#cData = getChain(single_ele_Run2015D, histname='')
cData = getChain([QCDHT_25ns, WJetsHTToLNu_25ns, TTJets_combined_2, singleTop_25ns, DY_25ns, TTV_25ns] , histname='')

#get template for fit method
template_QCD = ROOT.TH1F('template_QCD','template_QCD',30,-0.5,2.5)
templateName, templateCut = nameAndCut((250,-1), (500,-1), (3,4), (0,0), presel=antiSelStr, charge="", btagVar = 'nBJetMediumCSV30')
#cData.Draw('Lp>>template_QCD','('+templateCut+trigger+filters+')','goff')
cData.Draw('Lp>>template_QCD','('+weight_str+')*('+templateCut+')','goff')

histos = {}
bins = {}
fitRes = {}
#perform the fits in QCD CR
for crNJet in sorted(fitCR):
  fitRes[crNJet] = {}
  for ltb in sorted(fitCR[crNJet]):
    fitRes[crNJet][ltb] = {} 
    for htb in sorted(fitCR[crNJet][ltb]):
      deltaPhiCut = 1.0
      antiSelname, antiSelCut = nameAndCut(ltb, htb, crNJet, btb=(0,0), presel=antiSelStr, charge="", btagVar = 'nBJetMediumCSV30')
      Selname, SelCut         = nameAndCut(ltb, htb, crNJet, btb=(0,0), presel=SelStr, charge="", btagVar = 'nBJetMediumCSV30')      
      
      histos['QCD']={}
      histos['EWK']={}
      histos['DATA']={}
      histos['QCD']['antiSelection']=ROOT.TH1F('QCD_antiSelection','QCD_antiSelection',30,-0.5,2.5)
      histos['QCD']['Selection']=ROOT.TH1F('QCD_Selection','QCD_Selection',30,-0.5,2.5)
      histos['EWK']['antiSelection']=ROOT.TH1F('EWK_antiSelection','EWK_antiSelection',30,-0.5,2.5)
      histos['EWK']['Selection']=ROOT.TH1F('EWK_Selection','EWK_Selection',30,-0.5,2.5)
      histos['DATA']['antiSelection']=ROOT.TH1F('DATA_antiSelection','DATA_antiSelection',30,-0.5,2.5)
      histos['DATA']['Selection']=ROOT.TH1F('DATA_Selection','DATA_Selection',30,-0.5,2.5)

      Canv = ROOT.TCanvas('Canv','Canv')
      #mergeCanv.SetLogy()
      leg = ROOT.TLegend(0.7,0.75,0.98,0.95)
      leg.SetFillColor(0)
      leg.SetBorderSize(1)
      leg.SetShadowColor(ROOT.kWhite)
      text = ROOT.TLatex()
      text.SetNDC()
      text.SetTextSize(0.045)
      text.SetTextAlign(11)

      cQCD.Draw('Lp>>QCD_antiSelection','('+weight_str+')*('+antiSelCut+')')
      cQCD.Draw('Lp>>QCD_Selection','('+weight_str+')*('+SelCut+')')
      cEWK.Draw('Lp>>EWK_antiSelection','('+weight_str+')*('+antiSelCut+')')
      cEWK.Draw('Lp>>EWK_Selection','('+weight_str+')*('+SelCut+')')
#      cData.Draw('Lp>>DATA_antiSelection','('+antiSelCut+trigger+filters+')')
      cData.Draw('Lp>>DATA_antiSelection','('+weight_str+')*('+antiSelCut+')')
#      cData.Draw('Lp>>DATA_Selection','('+SelCut+trigger+filters+')')
      cData.Draw('Lp>>DATA_Selection','('+weight_str+')*('+SelCut+')')

#      rCSanti = getRCS(cData, antiSelCut+trigger+filters, deltaPhiCut, useWeight = False, weight = weight_str)
#      rCSsel = getRCS(cData, SelCut+trigger+filters, deltaPhiCut, useWeight = False, weight = weight_str)
      rCSanti = getRCS(cData, antiSelCut, deltaPhiCut, useWeight = True, weight = weight_str)
      rCSsel = getRCS(cData, SelCut, deltaPhiCut, useWeight = True, weight = weight_str)

      for hist in [histos['DATA']['antiSelection'],histos['DATA']['Selection']]:
        hist.SetStats(0)
        hist.GetYaxis().SetTitle('# of Events')
        hist.GetXaxis().SetTitle('L_{p}')
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineStyle(1)
        hist.SetLineWidth(1)

      for hist in [histos['QCD']['antiSelection'],histos['QCD']['Selection'],histos['EWK']['antiSelection'],histos['EWK']['Selection']]:
        hist.SetStats(0)
        hist.GetYaxis().SetTitle('# of Events')
        hist.GetXaxis().SetTitle('L_{p}')
        hist.SetLineWidth(2)
        hist.SetMarkerStyle(1)

      nEWKSel_err = ROOT.Double()
      nEWKSel = histos['EWK']['Selection'].IntegralAndError(0,histos['EWK']['Selection'].GetNbinsX(),nEWKSel_err)
      nEWKAntiSel_err = ROOT.Double()
      nEWKAntiSel = histos['EWK']['antiSelection'].IntegralAndError(0,histos['EWK']['antiSelection'].GetNbinsX(),nEWKAntiSel_err)
      nQCDSel_err = ROOT.Double()
      nQCDSel =  histos['QCD']['Selection'].IntegralAndError(0,histos['QCD']['Selection'].GetNbinsX(),nQCDSel_err)
      nQCDAntiSel_err = ROOT.Double()
      nQCDAntiSel = histos['QCD']['antiSelection'].IntegralAndError(0,histos['QCD']['antiSelection'].GetNbinsX(),nQCDAntiSel_err)
      nDATASel_err = ROOT.Double()
      nDATASel = histos['DATA']['Selection'].IntegralAndError(0,histos['DATA']['Selection'].GetNbinsX(),nDATASel_err)
      nDATAAntiSel_err = ROOT.Double()
      nDATAAntiSel = histos['DATA']['antiSelection'].IntegralAndError(0,histos['DATA']['antiSelection'].GetNbinsX(),nDATAAntiSel_err)

      fitRes[crNJet][ltb][htb] = {'NDATASel':nDATASel, 'NDATASel_err':float(nDATASel_err),\
                                  'NDATAAntiSel':nDATAAntiSel, 'NDATAAntiSel_err':float(nDATAAntiSel_err),\
                                  'NEWKSelMC':nEWKSel, 'NEWKSelMC_err':float(nEWKSel_err),\
                                  'NEWKAntiSelMC':nEWKAntiSel, 'NEWKAntiSelMC_err':float(nEWKAntiSel_err),\
                                  'NQCDSelMC':nQCDSel, 'NQCDSelMC_err':float(nQCDSel_err),\
                                  'NQCDAntiSelMC':nQCDAntiSel, 'NQCDAntiSelMC_err':float(nQCDAntiSel_err),\
                                  'deltaPhiCut':deltaPhiCut, 'rCSselectedDATA':rCSsel, 'rCSantiSelectedDATA':rCSanti}

      Canv.cd()
      histos['QCD']['antiSelection'].SetLineColor(ROOT.kRed)
      histos['QCD']['antiSelection'].SetLineStyle(ROOT.kDashed)
      leg.AddEntry(histos['QCD']['antiSelection'],'QCD anti-selected','l')

      histos['QCD']['Selection'].SetLineColor(ROOT.kRed)
      leg.AddEntry(histos['QCD']['Selection'],'QCD selected','l')

      histos['EWK']['antiSelection'].SetLineColor(ROOT.kBlack)
      histos['EWK']['antiSelection'].SetLineStyle(ROOT.kDashed)
      leg.AddEntry(histos['EWK']['antiSelection'],'EWK anti-selected','l')

      histos['EWK']['Selection'].SetLineColor(ROOT.kBlack)
      leg.AddEntry(histos['EWK']['Selection'],'EWK selected','l')

      histos['DATA']['antiSelection'].SetMarkerStyle(24)
      histos['DATA']['Selection'].SetMarkerStyle(20)
      leg.AddEntry(histos['DATA']['antiSelection'],'Data anti-selected')
      leg.AddEntry(histos['DATA']['Selection'],'Data selected')

      histos['QCD']['antiSelection'].Draw('hist e')
      histos['QCD']['Selection'].Draw('hist same e')
      histos['EWK']['antiSelection'].Draw('hist same e')
      histos['EWK']['Selection'].Draw('hist same e')
      histos['DATA']['antiSelection'].Draw('same ep')
      histos['DATA']['Selection'].Draw('same ep')

      histos['QCD']['antiSelection'].SetMaximum(1.5*histos['QCD']['antiSelection'].GetMaximum())
      histos['QCD']['Selection'].SetMaximum(1.5*histos['QCD']['Selection'].GetMaximum())
      histos['EWK']['antiSelection'].SetMaximum(1.5*histos['EWK']['antiSelection'].GetMaximum())
      histos['EWK']['Selection'].SetMaximum(1.5*histos['EWK']['Selection'].GetMaximum())
      histos['DATA']['antiSelection'].SetMaximum(1.5*histos['DATA']['antiSelection'].GetMaximum())
      histos['DATA']['Selection'].SetMaximum(1.5*histos['DATA']['Selection'].GetMaximum())

      leg.Draw()
#      text.DrawLatex(0.16,.96,"CMS #bf{#it{Preliminary}}")
      text.DrawLatex(0.16,.96,"CMS #bf{#it{Simulation}}")
      text.DrawLatex(0.62,0.96,"#bf{L="+str(lumi)+" fb^{-1} (13 TeV)}")

      Canv.cd()
      Canv.Print(wwwDir+prefix+Selname+'.png')
      Canv.Print(wwwDir+prefix+Selname+'.root')
      Canv.Print(wwwDir+prefix+Selname+'.pdf')
      Canv.Clear()

      LpTemplates = {'DATAantiSel':template_QCD, 'DATAsel':histos['DATA']['Selection'],\
                     'EWKantiSel':histos['EWK']['antiSelection'], 'EWKsel':histos['EWK']['Selection'],\
                     'QCDantiSel':histos['QCD']['antiSelection'], 'QCDsel':histos['QCD']['Selection']}
      fit_QCD = LpTemplateFit(LpTemplates, prefix=prefix+Selname, printDir=wwwDir+'templateFit')
      fitRes[crNJet][ltb][htb].update(fit_QCD)
      try: F_ratio = fit_QCD['QCD']['yield']/nDATAAntiSel
      except ZeroDivisionError: F_ratio = float('nan')
      try: F_ratio_err = F_ratio*sqrt(fit_QCD['QCD']['yieldVar']/fit_QCD['QCD']['yield']**2 + nDATAAntiSel_err**2/nDATAAntiSel**2)
      except ZeroDivisionError: F_ratio_err = float('nan')
      fitRes[crNJet][ltb][htb].update({'F_seltoantisel':F_ratio, 'F_seltoantisel_err':F_ratio_err})

      ROOT.setTDRStyle()
      if not os.path.exists(picklePath):
        os.makedirs(picklePath)
      pickle.dump(fitRes, file(picklePath+pickleFit,'w'))

for srNJet in sorted(signalRegion):
  bins[srNJet] = {}
  for stb in sorted(signalRegion[srNJet]):
    bins[srNJet][stb] = {}
    for htb in sorted(signalRegion[srNJet][stb]):
      bins[srNJet][stb][htb] = {}
      for btb in btreg:
        bins[srNJet][stb][htb][btb] = {}
        for dP in sorted(signalRegion[srNJet][stb][htb]):
          deltaPhiCut = signalRegion[srNJet][stb][htb][dP]['deltaPhi']

          print 'Binning => Ht: ',htb,'Lt: ',stb,'NJet: ',srNJet
          antiSelname, antiSelCut = nameAndCut(stb, htb, srNJet, btb=btb, presel=antiSelStr, charge="", btagVar = 'nBJetMediumCSV30')
          Selname, SelCut         = nameAndCut(stb, htb, srNJet, btb=btb, presel=SelStr, charge="", btagVar = 'nBJetMediumCSV30')

          histos['QCD']={}
          histos['EWK']={}
          histos['DATA']={}
          histos['QCD']['antiSelection']=ROOT.TH1F('QCD_antiSelection','QCD_antiSelection',30,-0.5,2.5)
          histos['QCD']['Selection']=ROOT.TH1F('QCD_Selection','QCD_Selection',30,-0.5,2.5)
          histos['EWK']['antiSelection']=ROOT.TH1F('EWK_antiSelection','EWK_antiSelection',30,-0.5,2.5)
          histos['EWK']['Selection']=ROOT.TH1F('EWK_Selection','EWK_Selection',30,-0.5,2.5)
          histos['DATA']['antiSelection']=ROOT.TH1F('DATA_antiSelection','DATA_antiSelection',30,-0.5,2.5)
          histos['DATA']['Selection']=ROOT.TH1F('DATA_Selection','DATA_Selection',30,-0.5,2.5)

          Canv = ROOT.TCanvas('Canv','Canv')
          #mergeCanv.SetLogy()
          leg = ROOT.TLegend(0.7,0.75,0.98,0.95)
          leg.SetFillColor(0)
          leg.SetBorderSize(1)
          leg.SetShadowColor(ROOT.kWhite)
          text = ROOT.TLatex()
          text.SetNDC()
          text.SetTextSize(0.045)
          text.SetTextAlign(11)
        
          cQCD.Draw('Lp>>QCD_antiSelection','('+weight_str+')*('+antiSelCut+')')
          cQCD.Draw('Lp>>QCD_Selection','('+weight_str+')*('+SelCut+')')
          cEWK.Draw('Lp>>EWK_antiSelection','('+weight_str+')*('+antiSelCut+')')
          cEWK.Draw('Lp>>EWK_Selection','('+weight_str+')*('+SelCut+')')
#          cData.Draw('Lp>>DATA_antiSelection','('+antiSelCut+trigger+filters+')')
          cData.Draw('Lp>>DATA_antiSelection','('+weight_str+')*('+antiSelCut+')')
#          cData.Draw('Lp>>DATA_Selection','('+SelCut+trigger+filters+')')
          cData.Draw('Lp>>DATA_Selection','('+weight_str+')*('+SelCut+')')

#          rCSanti = getRCS(cData, antiSelCut+trigger+filters, deltaPhiCut, useWeight = False, weight = weight_str)
#          rCSsel = getRCS(cData, SelCut+trigger+filters, deltaPhiCut, useWeight = False, weight = weight_str)
          rCSanti = getRCS(cData, antiSelCut, deltaPhiCut, useWeight = True, weight = weight_str)
          rCSsel = getRCS(cData, SelCut, deltaPhiCut, useWeight = True, weight = weight_str)

          for hist in [histos['DATA']['antiSelection'],histos['DATA']['Selection']]:
            hist.SetStats(0)
            hist.GetYaxis().SetTitle('# of Events')
            hist.GetXaxis().SetTitle('L_{p}')
            hist.SetLineColor(ROOT.kBlack)
            hist.SetLineStyle(1)
            hist.SetLineWidth(1)
  
          for hist in [histos['QCD']['antiSelection'],histos['QCD']['Selection'],histos['EWK']['antiSelection'],histos['EWK']['Selection']]:
            hist.SetStats(0)
            hist.GetYaxis().SetTitle('# of Events')
            hist.GetXaxis().SetTitle('L_{p}')
            hist.SetLineWidth(2)
            hist.SetMarkerStyle(1)

          nEWKSel_err = ROOT.Double()
          nEWKSel = histos['EWK']['Selection'].IntegralAndError(0,histos['EWK']['Selection'].GetNbinsX(),nEWKSel_err)
          nEWKAntiSel_err = ROOT.Double()
          nEWKAntiSel = histos['EWK']['antiSelection'].IntegralAndError(0,histos['EWK']['antiSelection'].GetNbinsX(),nEWKAntiSel_err)
          nQCDSel_err = ROOT.Double()
          nQCDSel =  histos['QCD']['Selection'].IntegralAndError(0,histos['QCD']['Selection'].GetNbinsX(),nQCDSel_err) 
          nQCDAntiSel_err = ROOT.Double()
          nQCDAntiSel = histos['QCD']['antiSelection'].IntegralAndError(0,histos['QCD']['antiSelection'].GetNbinsX(),nQCDAntiSel_err)
          nDATASel_err = ROOT.Double()
          nDATASel = histos['DATA']['Selection'].IntegralAndError(0,histos['DATA']['Selection'].GetNbinsX(),nDATASel_err)
          nDATAAntiSel_err = ROOT.Double()
          nDATAAntiSel = histos['DATA']['antiSelection'].IntegralAndError(0,histos['DATA']['antiSelection'].GetNbinsX(),nDATAAntiSel_err)

          bins[srNJet][stb][htb][btb][dP] = {'NDATASel':nDATASel, 'NDATASel_err':float(nDATASel_err),\
                                             'NDATAAntiSel':nDATAAntiSel, 'NDATAAntiSel_err':float(nDATAAntiSel_err),\
                                             'NEWKSelMC':nEWKSel, 'NEWKSelMC_err':float(nEWKSel_err),\
                                             'NEWKAntiSelMC':nEWKAntiSel, 'NEWKAntiSelMC_err':float(nEWKAntiSel_err),\
                                             'NQCDSelMC':nQCDSel, 'NQCDSelMC_err':float(nQCDSel_err),\
                                             'NQCDAntiSelMC':nQCDAntiSel, 'NQCDAntiSelMC_err':float(nQCDAntiSel_err),\
                                             'deltaPhiCut':deltaPhiCut, 'rCSselectedDATA':rCSsel, 'rCSantiSelectedDATA':rCSanti}
  
          Canv.cd()
          histos['QCD']['antiSelection'].SetLineColor(ROOT.kRed)
          histos['QCD']['antiSelection'].SetLineStyle(ROOT.kDashed)
          leg.AddEntry(histos['QCD']['antiSelection'],'QCD anti-selected','l')
   
          histos['QCD']['Selection'].SetLineColor(ROOT.kRed)
          leg.AddEntry(histos['QCD']['Selection'],'QCD selected','l')
   
          histos['EWK']['antiSelection'].SetLineColor(ROOT.kBlack)
          histos['EWK']['antiSelection'].SetLineStyle(ROOT.kDashed)
          leg.AddEntry(histos['EWK']['antiSelection'],'EWK anti-selected','l')
   
          histos['EWK']['Selection'].SetLineColor(ROOT.kBlack)
          leg.AddEntry(histos['EWK']['Selection'],'EWK selected','l')
  
          histos['DATA']['antiSelection'].SetMarkerStyle(24)
          histos['DATA']['Selection'].SetMarkerStyle(20)
          leg.AddEntry(histos['DATA']['antiSelection'],'Data anti-selected')
          leg.AddEntry(histos['DATA']['Selection'],'Data selected') 
        
          histos['QCD']['antiSelection'].Draw('hist e')
          histos['QCD']['Selection'].Draw('hist same e')
          histos['EWK']['antiSelection'].Draw('hist same e')
          histos['EWK']['Selection'].Draw('hist same e')
          histos['DATA']['antiSelection'].Draw('same ep')
          histos['DATA']['Selection'].Draw('same ep')
  
          histos['QCD']['antiSelection'].SetMaximum(1.5*histos['QCD']['antiSelection'].GetMaximum())
          histos['QCD']['Selection'].SetMaximum(1.5*histos['QCD']['Selection'].GetMaximum())
          histos['EWK']['antiSelection'].SetMaximum(1.5*histos['EWK']['antiSelection'].GetMaximum())
          histos['EWK']['Selection'].SetMaximum(1.5*histos['EWK']['Selection'].GetMaximum())
          histos['DATA']['antiSelection'].SetMaximum(1.5*histos['DATA']['antiSelection'].GetMaximum())
          histos['DATA']['Selection'].SetMaximum(1.5*histos['DATA']['Selection'].GetMaximum())
            
          leg.Draw()
#          text.DrawLatex(0.16,.96,"CMS #bf{#it{Preliminary}}")
          text.DrawLatex(0.16,.96,"CMS #bf{#it{Simulation}}")
          text.DrawLatex(0.62,0.96,"#bf{L="+str(lumi)+" fb^{-1} (13 TeV)}")
  
          Canv.cd()
          Canv.Print(wwwDir+prefix+Selname+'.png')
          Canv.Print(wwwDir+prefix+Selname+'.root')
          Canv.Print(wwwDir+prefix+Selname+'.pdf')
          Canv.Clear()
  
          ROOT.setTDRStyle()
          if not os.path.exists(picklePath):
            os.makedirs(picklePath)
          pickle.dump(bins, file(picklePath+picklePresel,'w'))

#derive N_QCD(dPhi<x) and N_QCD(dPhi>x)
for srNJet in sorted(signalRegion):
  for stb in sorted(signalRegion[srNJet]):
    for htb in sorted(signalRegion[srNJet][stb]):
      for btb in btreg:
        for dP in sorted(signalRegion[srNJet][stb][htb]):
          deltaPhiCut = signalRegion[srNJet][stb][htb][dP]['deltaPhi']
          sys         = signalRegion[srNJet][stb][htb][dP]['sys']
          Fsta        = fitRes[(3,4)][stb][(500,-1)]['F_seltoantisel']
          Fsta_err    = fitRes[(3,4)][stb][(500,-1)]['F_seltoantisel_err']
          Nanti       = bins[srNJet][stb][htb][btb][dP]['NDATAAntiSel']
          Nanti_err   = bins[srNJet][stb][htb][btb][dP]['NDATAAntiSel_err']
          RcsAnti     = bins[srNJet][stb][htb][btb][dP]['rCSantiSelectedDATA']['rCS']
          RcsAnti_err = bins[srNJet][stb][htb][btb][dP]['rCSantiSelectedDATA']['rCSE_pred']
          NQCD        = Fsta * Nanti
          NQCD_err    = sqrt( (Fsta_err**2*Nanti**2+Nanti_err**2*Fsta**2) + (sys)**2 )
          try: NQCD_lowDPhi = NQCD/(RcsAnti+1)
          except ZeroDivisionError: NQCD_lowDPhi = float('nan') 
          try: NQCD_lowDPhi_err = NQCD_lowDPhi*sqrt((NQCD_err/NQCD)**2 + (RcsAnti_err/(RcsAnti+1))**2)
          except ZeroDivisionError: NQCD_lowDPhi_err = float('nan')
          NQCD_highDPhi = NQCD - NQCD_lowDPhi
          NQCD_highDPhi_err = sqrt(RcsAnti_err**2*NQCD_lowDPhi**2 + RcsAnti**2*NQCD_lowDPhi_err**2)
          bins[srNJet][stb][htb][btb][dP].update({'NQCDpred':NQCD, 'NQCDpred_err':NQCD_err, 'NQCDpred_lowdPhi':NQCD_lowDPhi, 'NQCDpred_lowdPhi_err':NQCD_lowDPhi_err, 'NQCDpred_highdPhi':NQCD_highDPhi, 'NQCDpred_highdPhi_err':NQCD_highDPhi_err})
          pickle.dump(bins, file(picklePath+picklePresel,'w'))

