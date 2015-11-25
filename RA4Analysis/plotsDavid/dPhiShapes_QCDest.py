import ROOT
import pickle 
import copy, os, sys
ROOT.gROOT.LoadMacro("../../HEPHYPythonTools/scripts/root/tdrstyle.C")
ROOT.TH1F().SetDefaultSumw2()
ROOT.setTDRStyle()
ROOT.gStyle.SetMarkerStyle(20)
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

preprefix = 'QCDestimation/newAntiID/hOverE/normalization/'
wwwDir = '/afs/hephy.at/user/'+username[0]+'/'+username+'/www/RunII/Spring15_25ns/'+preprefix+'/'
prefix = 'Lp_singleElectronic'
 
if not os.path.exists(wwwDir):
  os.makedirs(wwwDir)

signalRegion = {(3, 4): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}, #3-4jets QCD and W+jets control region
                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (350, -1):  {(500, -1):   {(0.75):   {'deltaPhi': 0.75}}},
                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0},
                                                    (0.75):   {'deltaPhi': 0.75}},
                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0},
                                                    (0.75):   {'deltaPhi': 0.75}},
                                      (500, 750):  {(0.75):   {'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'deltaPhi': 0.75}}}},
                (4, 5): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}, #4-5jets TTbar control region
                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (350, -1):  {(500, -1):   {(0.75):   {'deltaPhi': 0.75}}},
                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0},
                                                    (0.75):   {'deltaPhi': 0.75}},
                                      (500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0},
                                                    (0.75):   {'deltaPhi': 0.75}},
                                      (500, 750):  {(0.75):   {'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'deltaPhi': 0.75}}}},
                (5, 5): {(250, 350): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},  #signal regions
                         (350, 450): {(500, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (450, -1):  {(500, -1):   {(1.0):    {'deltaPhi': 1.0},
                                                    (0.75):   {'deltaPhi': 0.75}}}},
                (6, 7): {(250, 350): {(500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                         (350, 450): {(500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                          (450, -1): {(500, 750):  {(0.75):   {'deltaPhi': 0.75}},
                                      (750, -1):   {(0.75):   {'deltaPhi': 0.75}},
                                      (500, 1000): {(0.75):   {'deltaPhi': 0.75}},
                                      (1000, -1):  {(0.75):   {'deltaPhi': 0.75}}}},
                (8, -1): {(250, 350):{(500, 750):  {(1.0):    {'deltaPhi': 1.0}},
                                      (750, -1):   {(1.0):    {'deltaPhi': 1.0}}},
                          (350, -1): {(500, -1):   {(0.75):   {'deltaPhi': 0.75}}},
                          (350, 450):{(500, -1):   {(0.75):   {'deltaPhi': 0.75}}},
                          (450, -1): {(500, -1):   {(0.75):   {'deltaPhi': 0.75}}}}
}

btreg = [(0,0)]#, (1,1), (2,-1)] #1b and 2b estimates are needed for the btag fit

def makeWeight(lumi=4., sampleLumi=3.,debug=False):
  if debug:
    print 'No lumi-reweighting done!!'
    return 'weight', 'weight*weight'
  else:
    weight_str = '(((weight)/'+str(sampleLumi)+')*'+str(lumi)+')'
    weight_err_str = '('+weight_str+'*'+weight_str+')'
    return weight_str, weight_err_str
lumi = 3.
sampleLumi = 1.55
debugReweighting = False
weight_str, weight_err_str = makeWeight(lumi, sampleLumi=sampleLumi, debug=debugReweighting)

def getRCS(c, cut, dPhiCut, useWeight = False, weight = 'weight'):
#  dPhiStr = 'acos((LepGood_pt+met_pt*cos(LepGood_phi-met_phi))/sqrt(LepGood_pt**2+met_pt**2+2*met_pt*LepGood_pt*cos(LepGood_phi-met_phi)))'
  dPhiStr = 'deltaPhi_Wl'
  if useWeight:
    h = getPlotFromChain(c, dPhiStr, [0,dPhiCut,pi], cutString=cut, binningIsExplicit=True, weight = weight)
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

#trigger and filters for real Data
trigger = "&&(HLT_EleHT350||HLT_MuHT350)"
filters = "&&Flag_CSCTightHaloFilter&&Flag_HBHENoiseFilter_fix&&Flag_HBHENoiseFilter&&Flag_goodVertices&&Flag_eeBadScFilter&&Flag_EcalDeadCellTriggerPrimitiveFilter"
#filters = "&&Flag_CSCTightHaloFilter&&Flag_HBHENoiseFilter_fix&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_eeBadScFilter"

presel = 'nLep==1&&nVeto==0&&nEl==1&&leptonPt>25&&Jet2_pt>80'
antiSelStr = presel+filters+'&&Selected==-1&&leptonHoverE>0.01'
SelStr = presel+filters+'&&Selected==1'

cQCD  = getChain(QCDHT_25ns,histname='')
#cEWK  = getChain([WJetsHTToLNu_25ns, TTJets_HTLO_25ns, singleTop_25ns, DY_25ns, TTV_25ns],histname='')
#cData = getChain(data_ele_25ns , histname='')

histos = {}
bins = {}
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
#          histos['EWK']={}
#          histos['DATA']={}
          histos['QCD']['antiSelection']=ROOT.TH1F('QCD_antiSelection','QCD_antiSelection',30,-0.5,2.5)
          histos['QCD']['Selection']=ROOT.TH1F('QCD_Selection','QCD_Selection',30,-0.5,2.5)
#          histos['EWK']['antiSelection']=ROOT.TH1F('EWK_antiSelection','EWK_antiSelection',30,-0.5,2.5)
#          histos['EWK']['Selection']=ROOT.TH1F('EWK_Selection','EWK_Selection',30,-0.5,2.5)
#          histos['DATA']['antiSelection']=ROOT.TH1F('DATA_antiSelection','DATA_antiSelection',30,-0.5,2.5)
#          histos['DATA']['Selection']=ROOT.TH1F('DATA_Selection','DATA_Selection',30,-0.5,2.5)

          Canv = ROOT.TCanvas('Canv','Canv')
#          Canv.SetLogy()
          leg = ROOT.TLegend(0.7,0.8,0.98,0.95)
          leg.SetFillColor(0)
          leg.SetBorderSize(1)
          leg.SetShadowColor(ROOT.kWhite)
          text = ROOT.TLatex()
          text.SetNDC()
          text.SetTextSize(0.045)
          text.SetTextAlign(11)
        
          cQCD.Draw('Lp>>QCD_antiSelection','('+weight_str+')*('+antiSelCut+')')
          cQCD.Draw('Lp>>QCD_Selection','('+weight_str+')*('+SelCut+')')
#          cEWK.Draw('Lp>>EWK_antiSelection','(weight)*('+antiSelCut+')')
#          cEWK.Draw('Lp>>EWK_Selection','(weight)*('+SelCut+')')
#          cData.Draw('Lp>>DATA_antiSelection','('+antiSelCut+trigger+')')
#          cData.Draw('Lp>>DATA_Selection','('+SelCut+trigger+')')

#          for hist in [histos['DATA']['antiSelection'],histos['DATA']['Selection']]:
#            hist.SetStats(0)
#            hist.GetYaxis().SetTitle('# of Events')
#            hist.GetXaxis().SetTitle('L_{p}')
#            hist.SetLineColor(ROOT.kBlack)
#            hist.SetLineStyle(1)
#            hist.SetLineWidth(1)
  
          for hist in [histos['QCD']['antiSelection'],histos['QCD']['Selection']]:#,histos['EWK']['antiSelection'],histos['EWK']['Selection']]:
            hist.SetStats(0)
            hist.GetYaxis().SetTitle('a.u.')
            hist.GetXaxis().SetTitle('L_{p}')
            hist.SetLineWidth(2)
            hist.SetMarkerStyle(20)

          Canv.cd()
          histos['QCD']['antiSelection'].SetLineColor(ROOT.kBlue)
          histos['QCD']['antiSelection'].SetMarkerColor(ROOT.kBlue)
#          histos['QCD']['antiSelection'].SetLineStyle(ROOT.kDashed)
          if histos['QCD']['antiSelection'].Integral()>0: histos['QCD']['antiSelection'].Scale(1./histos['QCD']['antiSelection'].Integral())
          leg.AddEntry(histos['QCD']['antiSelection'],'QCD anti-selected','l')
   
          histos['QCD']['Selection'].SetLineColor(ROOT.kRed)
          histos['QCD']['Selection'].SetMarkerColor(ROOT.kRed)
          if histos['QCD']['Selection'].Integral()>0: histos['QCD']['Selection'].Scale(1./histos['QCD']['Selection'].Integral())
          leg.AddEntry(histos['QCD']['Selection'],'QCD selected','l')
   
          histos['QCD']['antiSelection'].Draw('ep')
          histos['QCD']['Selection'].Draw('same ep')
  
          histos['QCD']['antiSelection'].SetMaximum(1.5*histos['QCD']['antiSelection'].GetMaximum())
          histos['QCD']['Selection'].SetMaximum(1.5*histos['QCD']['Selection'].GetMaximum())
#          histos['QCD']['antiSelection'].SetMinimum(0)
#          histos['QCD']['Selection'].SetMinimum(0)
            
          leg.Draw()
          text.DrawLatex(0.16,.96,"CMS #bf{#it{Simulation}}")
          text.DrawLatex(0.75,0.96,"#bf{MC (13 TeV)}")
  
          Canv.cd()
          Canv.Print(wwwDir+prefix+Selname+'.png')
          Canv.Print(wwwDir+prefix+Selname+'.root')
          Canv.Print(wwwDir+prefix+Selname+'.pdf')
          Canv.Clear()
  
