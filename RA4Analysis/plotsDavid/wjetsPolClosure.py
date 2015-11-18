import ROOT
import pickle 
import copy, os, sys
ROOT.gROOT.LoadMacro("../../HEPHYPythonTools/scripts/root/tdrstyle.C")
ROOT.TH1F().SetDefaultSumw2()
ROOT.setTDRStyle()
#ROOT.tdrStyle.SetPadBottomMargin(0.2)
#ROOT.tdrStyle.SetPadRightMargin(0.25)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetHistMinimumZero()

from Workspace.HEPHYPythonTools.helpers import *
from Workspace.RA4Analysis.helpers import *
from Workspace.RA4Analysis.signalRegions import *
from draw_helpers import *
from math import *
from Workspace.HEPHYPythonTools.user import username

preprefix = 'WPolarizationEstimation/closurePlots'
wwwDir = '/afs/hephy.at/user/d/dhandl/www/RunII/Spring15_25ns/'+preprefix+'/'
presel = ''

if not os.path.exists(wwwDir):
  os.makedirs(wwwDir)

path = '/data/'+username+'/Spring15/25ns/PredictionAN_3.0/'
pickleFileWPolPlus10  = 'singleLeptonic_Spring15_WPolPlus10_estimationResults_pkl_kappa_corrected'
pickleFileWPolMinus10 = 'singleLeptonic_Spring15_WPolMinus10_estimationResults_pkl_kappa_corrected'
pickleFileTTPolPlus5  = 'singleLeptonic_Spring15_TTPolPlus5_estimationResults_pkl_kappa_corrected'
pickleFileTTPolMinus5 = 'singleLeptonic_Spring15_TTPolMinus5_estimationResults_pkl_kappa_corrected'
pickleFile            = 'singleLeptonic_Spring15_estimationResults_pkl_kappa_corrected'
resWPolPlus10 = pickle.load(file(path+pickleFileWPolPlus10))
resWPolMinus10 = pickle.load(file(path+pickleFileWPolMinus10))
resTTPolPlus5 = pickle.load(file(path+pickleFileTTPolPlus5))
resTTPolMinus5 = pickle.load(file(path+pickleFileTTPolMinus5))
resTruth = pickle.load(file(path+pickleFile))

targetLumi = 3.

def getFraction(Bkg, Bkg_err, QCD, QCD_err):
  try: res = QCD/Bkg
  except ZeroDivisionError: res = float('nan')
  try: res_err = res*sqrt(Bkg_err**2/Bkg**2 + QCD_err**2/QCD**2)
  except ZeroDivisionError: res_err = float('nan')
  return res, res_err

#define SR
signalRegions = signalRegion3fb
signalRegion = {#(3, 4): {(250, 350): {(500, -1):   {'deltaPhi': 1.0}, #3-4jets QCD and W+jets control region
                #                      (500, 750):  {'deltaPhi': 1.0},
                #                      (750, -1):   {'deltaPhi': 1.0}},
                #         (350, 450): {(500, -1):   {'deltaPhi': 1.0},
                #                      (500, -1):   {'deltaPhi': 0.75},
                #                      (500, 750):  {'deltaPhi': 1.0},
                #                      (750, -1):   {'deltaPhi': 1.0}},
                #         (450, -1):  {(500, -1):   {'deltaPhi': 1.0},
                #                      (500, -1):   {'deltaPhi': 0.75},
                #                      (500, 1000): {'deltaPhi': 0.75},
                #                      (1000, -1):  {'deltaPhi': 0.75}}},
                #(4, 5): {(250, 350): {(500, -1):   {'deltaPhi': 1.0}, #4-5jets TTbar control region
                #                      (500, 750):  {'deltaPhi': 1.0},
                #                      (750, -1):   {'deltaPhi': 1.0}},
                #         (350, 450): {(500, -1):   {'deltaPhi': 1.0},
                #                      (500, -1):   {'deltaPhi': 0.75},
                #                      (500, 750):  {'deltaPhi': 1.0},
                #                      (750, -1):   {'deltaPhi': 1.0}},
                #         (450, -1):  {(500, -1):   {'deltaPhi': 1.0},
                #                      (500, -1):   {'deltaPhi': 0.75},
                #                      (500, 1000): {'deltaPhi': 0.75},
                #                      (1000, -1):  {'deltaPhi': 0.75}}},
                (5, 5): {(250, 350): {(500, -1):   {'deltaPhi': 1.0}},  #signal regions
                         (350, 450): {(500, -1):   {'deltaPhi': 1.0}},
                         (450, -1):  {(500, -1):   {'deltaPhi': 1.0}}},
                (6, 7): {(250, 350): {(500, 750):  {'deltaPhi': 1.0},
                                      (750, -1):   {'deltaPhi': 1.0}},
                         (350, 450): {(500, 750):  {'deltaPhi': 1.0},
                                      (750, -1):   {'deltaPhi': 1.0}},
                          (450, -1): {(500, 1000): {'deltaPhi': 0.75},
                                      (1000, -1):  {'deltaPhi': 0.75}}},
                (8, -1): {(250, 350):{(500, 750):  {'deltaPhi': 1.0},
                                      (750, -1):   {'deltaPhi': 1.0}},
                          (350, 450):{(500, -1):   {'deltaPhi': 0.75}},
                          (450, -1): {(500, -1):   {'deltaPhi': 0.75}}}
}

btreg = [(0,0)]#, (1,1), (2,2)] #1b and 2b estimates are needed for the btag fit

rowsNJet = {}
rowsSt = {}
bins = 0
for srNJet in sorted(signalRegions):
  rowsNJet[srNJet] = {}
  rowsSt[srNJet] = {}
  rows = 0
  for stb in sorted(signalRegions[srNJet]):
    rows += len(signalRegions[srNJet][stb])
    rowsSt[srNJet][stb] = {'n':len(signalRegions[srNJet][stb])}
  rowsNJet[srNJet] = {'nST':len(signalRegions[srNJet]), 'n':rows}
  bins += rows

WPolPlus10  = []
WPolMinus10 = []
TTPolPlus5  = []
TTPolMinus5 = []
for srNJet in sorted(signalRegions):
  for stb in sorted(signalRegions[srNJet]):
    for htb in sorted(signalRegions[srNJet][stb]):
      uWPolPlus10  = ((resWPolPlus10[srNJet][stb][htb]['tot_pred']/resWPolPlus10[srNJet][stb][htb]['tot_truth'])/(resTruth[srNJet][stb][htb]['tot_pred']/resTruth[srNJet][stb][htb]['tot_truth'])) - 1.
      uWPolMinus10 = ((resWPolMinus10[srNJet][stb][htb]['tot_pred']/resWPolMinus10[srNJet][stb][htb]['tot_truth'])/(resTruth[srNJet][stb][htb]['tot_pred']/resTruth[srNJet][stb][htb]['tot_truth'])) - 1.
      uTTPolPlus5  = ((resTTPolPlus5[srNJet][stb][htb]['tot_pred']/resTTPolPlus5[srNJet][stb][htb]['tot_truth'])/(resTruth[srNJet][stb][htb]['tot_pred']/resTruth[srNJet][stb][htb]['tot_truth'])) - 1.
      uTTPolMinus5 = ((resTTPolMinus5[srNJet][stb][htb]['tot_pred']/resTTPolMinus5[srNJet][stb][htb]['tot_truth'])/(resTruth[srNJet][stb][htb]['tot_pred']/resTruth[srNJet][stb][htb]['tot_truth'])) - 1.
      WPolPlus10.append(uWPolPlus10)
      WPolMinus10.append(uWPolMinus10)
      TTPolPlus5.append(uTTPolPlus5)
      TTPolMinus5.append(uTTPolMinus5)

WPolPlus10_H  = ROOT.TH1F('WPolPlus10_H','WPolPlus10_H',bins,0,bins)
WPolMinus10_H = ROOT.TH1F('WPolMinus10_H','WPolMinus10_H',bins,0,bins)
WPolPlus10_H.SetLineColor(0)
WPolMinus10_H.SetLineColor(0)
WPolPlus10_H.SetFillColor(color('wJets'))
WPolMinus10_H.SetFillColor(color('wJets'))
WPolPlus10_H.SetStats(0)
WPolMinus10_H.SetStats(0)
WPolPlus10_H.SetMinimum(-0.05)
WPolMinus10_H.SetMinimum(-0.05)
WPolPlus10_H.SetMaximum(0.05)
WPolMinus10_H.SetMaximum(0.05)

TTPolPlus5_H  = ROOT.TH1F('TTPolPlus5_H','TTPolPlus5_H',bins,0,bins)
TTPolMinus5_H = ROOT.TH1F('TTPolMinus5_H','TTPolMinus5_H',bins,0,bins)
TTPolPlus5_H.SetLineColor(0)
TTPolMinus5_H.SetLineColor(0)
TTPolPlus5_H.SetFillColor(color('ttJets'))
TTPolMinus5_H.SetFillColor(color('ttJets'))
TTPolPlus5_H.SetStats(0)
TTPolMinus5_H.SetStats(0)
TTPolPlus5_H.SetMinimum(-0.05)
TTPolMinus5_H.SetMinimum(-0.05)
TTPolPlus5_H.SetMaximum(0.05)
TTPolMinus5_H.SetMaximum(0.05)

#w_pred_H  = ROOT.TH1F('w_pred_H','W+Jets pred.', bins,0,bins)
#w_truth_H = ROOT.TH1F('w_truth_H','W+Jets truth', bins,0,bins)
#w_pred_H.SetLineColor(color('wJets')+1)
#w_pred_H.SetFillColor(color('wJets'))
#w_pred_H.SetLineWidth(2)
#w_truth_H.SetLineColor(color('wJets')-1)
#w_truth_H.SetLineWidth(2)
#
#rest_H = ROOT.TH1F('rest_H','EWK rest', bins,0,bins)
#rest_H.SetLineColor(color('TTVH')+1)
#rest_H.SetFillColor(color('TTVH'))
#w_truth_H.SetLineColor(color('DY'))
#
#pred_H  = ROOT.TH1F('pred_H','total pred.', bins,0,bins)
#pred_H.SetBarWidth(0.4)
#pred_H.SetBarOffset(0.1)
#truth_H = ROOT.TH1F('truth_H','Total MC truth',bins,0,bins)
#truth_H.SetBarWidth(0.4)
#truth_H.SetBarOffset(0.1)
#pred_H.SetLineColor(ROOT.kGray+1)
#pred_H.SetMarkerStyle(1)
#pred_H.SetLineWidth(2)
#truth_H.SetLineColor(ROOT.kBlack)
#truth_H.SetLineWidth(2)
#
#drawOption = 'hist ][ e1'
#drawOptionSame = drawOption + 'same'
#
#predXErr = []
#predYErr = []
#predX = []
#predY = []

#i=1
#for srNJet in sorted(signalRegions):
#  for stb in sorted(signalRegions[srNJet]):
#    for htb in sorted(signalRegions[srNJet][stb]):
#      #print 1
#      #tt_pred_H.SetBinContent(i, resTruth[srNJet][stb][htb]['TT_pred'])
#      #tt_pred_H.SetBinError(i,   resTruth[srNJet][stb][htb]['TT_pred_err'])
#      #tt_truth_H.SetBinContent(i,resTruth[srNJet][stb][htb]['TT_truth'])
#      #tt_truth_H.SetBinError(i,  resTruth[srNJet][stb][htb]['TT_truth_err'])
#      
#      w_pred_H.SetBinContent(i, resWpolSys[srNJet][stb][htb]['W_pred'])
#      w_pred_H.SetBinError(i,   resWpolSys[srNJet][stb][htb]['W_pred_err'])
#      w_truth_H.SetBinContent(i,resWpolSysSR[srNJet][stb][htb]['W_truth'])
#      w_truth_H.SetBinError(i,  resWpolSysSR[srNJet][stb][htb]['W_truth_err'])
#      u = (resWpolSys[srNJet][stb][htb]['tot_pred']/(resWpolSys[srNJet][stb][htb]['tot_truth']+resWpolSysSR[srNJet][stb][htb]['W_truth']-resWpolSys[srNJet][stb][htb]['W_truth']))/(resTruth[srNJet][stb][htb]['tot_pred']/resTruth[srNJet][stb][htb]['tot_truth']) - 1.
#      predYErr.append(u)
#      predXErr.append(0.5)
#      #predY.append(resWpolSys[srNJet][stb][htb]['W_pred'])
#      predY.append(0.)
#      predX.append(i-0.5)
#
#       
#      #rest_H.SetBinContent(i,resWpolSys[srNJet][stb][htb]['Rest_truth'])
#      #rest_H.SetBinError(i,  resWpolSys[srNJet][stb][htb]['Rest_truth_err'])
#
#      #pred_H.SetBinContent(i, res[srNJet][stb][htb]['tot_pred'])
#      #pred_H.SetBinError(i,   res[srNJet][stb][htb]['tot_pred_err'])
#      #predYErr.append(res[srNJet][stb][htb]['tot_pred_err'])
#      #predXErr.append(0.5)
#      #predY.append(res[srNJet][stb][htb]['tot_pred'])
#      #predX.append(i-0.5)
#      #truth_H.SetBinContent(i,res[srNJet][stb][htb]['tot_truth'])
#      #truth_H.SetBinError(i,  res[srNJet][stb][htb]['tot_truth_err'])
#      #truth_H.GetXaxis().SetBinLabel(i, str(i))
#      #pred_H.GetXaxis().SetBinLabel(i, str(i))
#      i+=1
#
#ax = array('d',predX)
#ay = array('d',predY)
#aexh = array('d',predXErr)
#aexl = array('d',predXErr)
#aeyh = array('d',predYErr)
#aeyl = array('d',predYErr)

#can = ROOT.TCanvas('can','can',600,600)

#pad1=ROOT.TPad("pad1","MyTitle",0.,0.3,1.,1.)
#pad1=ROOT.TPad("pad1","MyTitle",0.,0.,1.,1.)
#pad1.SetLeftMargin(0.15)
#pad1.SetBottomMargin(0.)
#pad1.Draw()
#pad1.cd()

#h_Stack = ROOT.THStack('h_Stack','Stack')
#h_Stack.Add(tt_pred_H)
#h_Stack.Add(w_pred_H)
#h_Stack.Add(rest_H)
#h_Stack.SetMaximum(30)
#h_Stack.GetYaxis().SetTitle('Signal Region #')

#w_truth_H.SetMaximum(20)
#w_truth_H.GetYaxis().SetTitle('Events')
#w_truth_H.GetXaxis().SetTitle('Signal Region #')
#w_truth_H.GetYaxis().SetTitleSize(0.06)
#w_truth_H.GetYaxis().SetLabelSize(0.06)

#truth_H.SetMaximum(20)
#truth_H.GetYaxis().SetTitle('Events')
#truth_H.GetXaxis().SetTitle('Signal Region #')
#truth_H.GetYaxis().SetTitleSize(0.06)
#truth_H.GetYaxis().SetLabelSize(0.06)

#leg = ROOT.TLegend(0.65,0.7,0.98,0.95)
#leg.SetFillColor(ROOT.kWhite)
#leg.SetShadowColor(ROOT.kWhite)
#leg.SetBorderSize(1)
#leg.SetTextSize(0.045)
#leg.AddEntry(w_truth_H)
#leg.AddEntry(truth_H)
#leg.AddEntry(tt_pred_H,'','f')
#leg.AddEntry(w_pred_H,'','f')
#leg.AddEntry(rest_H,'','f')

#h_Stack.Draw('hist')
#h_Stack.GetYaxis().SetTitle('# of Events')
#h_Stack.GetYaxis().SetTitleOffset(0.8)
#h_Stack.GetYaxis().SetNdivisions(508)
#predError = ROOT.TGraphError(pred_H)
#pred_H.Draw('e1 same')
#pred_err = ROOT.TGraphAsymmErrors(bins, ax, ay, aexl, aexh, aeyl, aeyh)
#pred_err.SetFillColor(ROOT.kGray+1)
#pred_err.SetFillStyle(3244)
#pred_err.Draw('2 same')
#truth_H.Draw('e1p same')

#leg.Draw()

#can.cd()

#ratio = ROOT.TH1F('ratio','ratio',bins,0,bins)
#ratio.Sumw2()
#ratio = pred_H.Clone()
#ratio.Divide(truth_H)
#
#pad2=ROOT.TPad("pad2","datavsMC",0.,0.,1.,.3)
#pad2.SetLeftMargin(0.15)
#pad2.SetBottomMargin(0.3)
#pad2.SetTopMargin(0.)
#pad2.SetGrid()
#pad2.Draw()
#pad2.cd()
#ratio.SetLineColor(ROOT.kBlack)
#ratio.SetMarkerStyle(8)
#ratio.GetXaxis().SetTitle('Signal Region #')
#ratio.GetXaxis().SetTitleSize(0.13)
#ratio.GetXaxis().SetLabelSize(0.21)
#ratio.GetXaxis().SetNdivisions(508)
#ratio.GetYaxis().SetTitle('pred./truth')
#ratio.GetYaxis().SetTitleSize(0.13)
#ratio.GetYaxis().SetLabelSize(0.13)
#ratio.GetYaxis().SetTitleOffset(0.4)
#ratio.GetYaxis().SetNdivisions(508)
#ratio.SetMinimum(0.)
#ratio.SetMaximum(2.2)
#ratio.Draw('e1p')

#can.cd()

#ROOT_colors = [ROOT.kBlack, ROOT.kRed-4, ROOT.kBlue, ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kAzure+6, ROOT.kCyan+3, ROOT.kOrange , ROOT.kRed-10]
text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.04)
text.SetTextAlign(11)
canv = ROOT.TCanvas('canv','canv',600,600)
canv.SetGrid()
l = ROOT.TLegend(0.65,0.8,0.98,0.95)
l.SetFillColor(0)
l.SetBorderSize(1)
l.SetShadowColor(ROOT.kWhite)


j=0
for i_njb, njb in enumerate(sorted(signalRegion)):
  for i_CR, ltb in enumerate(sorted(signalRegion[njb])):
    for i_htb,htb in enumerate(sorted(signalRegion[njb][ltb])):
      j+=1
      WPolPlus10_H.SetBinContent(j,WPolPlus10[j-1])
      WPolMinus10_H.SetBinContent(j,WPolMinus10[j-1])
      TTPolPlus5_H.SetBinContent(j,TTPolPlus5[j-1])
      TTPolMinus5_H.SetBinContent(j,TTPolMinus5[j-1])
      TTPolMinus5_H.GetXaxis().SetBinLabel(j,str(j))

WPolPlus10_H.SetBarOffset(0.)
WPolPlus10_H.SetBarWidth(0.5)
WPolPlus10_H.Draw("bar")
WPolMinus10_H.SetBarOffset(0.)
WPolMinus10_H.SetBarWidth(0.5)
WPolMinus10_H.Draw("bar same")
TTPolPlus5_H.SetBarOffset(0.5)
TTPolPlus5_H.SetBarWidth(0.5)
TTPolPlus5_H.Draw("'bar same")
TTPolMinus5_H.SetBarOffset(0.5)
TTPolMinus5_H.SetBarWidth(0.5)
TTPolMinus5_H.GetXaxis().SetTitle('signal region')
TTPolMinus5_H.Draw("bar same")
l.AddEntry(WPolPlus10_H,'W+jets #pm10% var.','f')
l.AddEntry(TTPolPlus5_H,'t#bar{t}+jets #pm5% var.','f')
text.DrawLatex(0.17,.96,"CMS #bf{#it{Simulation}}")
text.DrawLatex(0.7,0.96,"#bf{L="+str(targetLumi)+" fb^{-1} (13 TeV)}")
l.Draw()
canv.Print(wwwDir+'WpolSystematics_totalBkg_inSR.png')
canv.Print(wwwDir+'WpolSystematics_totalBkg_inSR.pdf')
canv.Print(wwwDir+'WpolSystematics_totalBkg_inSR.root')
#
#canv2 = ROOT.TCanvas('canv2','canv2',600,600)
#
#ClosureHist=ROOT.TH1F('ClosureHist','ClosureHist',13,0,13)
#ClosureHist.SetLineWidth(2)
#k=0
#for i_njb, njb in enumerate(sorted(signalRegion)):
#  for i_CR, ltb in enumerate(sorted(signalRegion[njb])):
#    for i_htb,htb in enumerate(sorted(signalRegion[njb][ltb])):
#      k+=1
#      #res, res_err = getFraction(bins[(3,4)][ltb][htb][(0,0)]['NQCDSelMC'],bins[(3,4)][ltb][htb][(0,0)]['NQCDSelMC_err'],bins[(3,4)][ltb][htb][(0,0)]['NQCDpred'],bins[(3,4)][ltb][htb][(0,0)]['NQCDpred_err'])
#      res = (binsWpolSys[njb][ltb][htb]['W_pred']-bins[njb][ltb][htb]['W_pred'])/binsWpolSys[njb][ltb][htb]['W_pred']
#      ClosureHist.SetBinContent(k,res)
#      #ClosureHist.SetBinError(k,res_err)
#      #ClosureHist.GetXaxis().SetLabelSize(0.035)
#      #ClosureHist.GetXaxis().SetBinLabel(k,'nJ'+str(njb)+'_LT'+str(ltb)+'_HT'+str(htb))
#      ClosureHist.GetXaxis().SetBinLabel(k,k)
#      #ClosureHist.GetYaxis().SetTitle('#frac{N^{pred}_{}}{N^{MC}_{QCD}}')
#
#ClosureHist.Draw('P')
#ClosureHist.SetMinimum(-0.25)
#ClosureHist.SetMaximum(0.25)
#text.DrawLatex(0.15,.96,"CMS #bf{#it{Simulation}}")
#text.DrawLatex(0.6,0.96,"#bf{L="+str(targetLumi)+" fb^{-1} (13 TeV)}")
#line2 = ROOT.TLine()
#line2.SetY1(0.)
#line2.SetX2(13)
#line2.SetHorizontal()
#line2.SetLineColor(ROOT.kBlack)
#line2.SetLineStyle(ROOT.kDashed)
#line2.Draw()
#canv2.Print(wwwDir+presel+'WpolClosure_Wjets_inSR.png')
#canv2.Print(wwwDir+presel+'WpolClosure_Wjets_inSR.pdf')
#canv2.Print(wwwDir+presel+'WpolClosure_Wjets_inSR.root')

#plot F_sel-to-antisel binned in HT for all Njets
#ratio_ht={}
#for stb in streg:
#  ratio_ht[stb]={}
#  first = True
#  canv = ROOT.TCanvas('canv','canv',600,600)
#  #canv.SetLogy()
#  l = ROOT.TLegend(0.65,0.85,0.95,0.95)
#  l.SetFillColor(0)
#  l.SetBorderSize(1)
#  l.SetShadowColor(ROOT.kWhite)
#  
#  t=ROOT.TLatex()
#  t.SetNDC()
#  t.SetTextSize(0.04)
#  t.SetTextAlign(11)
#  for i_njb, njb in enumerate(njreg):
#    ratio_ht[stb][njb]={}
#    for btb in btreg:
#      ratio_ht[stb][njb][btb]=ROOT.TH1F('ratio_htHist','ratio_htHist',len(htreg),0,len(htreg))
#      ratio_ht[stb][njb][btb].SetLineColor(ROOT_colors[i_njb])
#      ratio_ht[stb][njb][btb].SetLineWidth(2)
#      for i_htb, htb in enumerate(htreg):
#        nQCDsel = bins[njb][stb][htb][btb]['nQCDSelected'] 
#        nQCDselVar = bins[njb][stb][htb][btb]['nQCDSelectedVar'] 
#        nQCDantisel = bins[njb][stb][htb][btb]['nAntiSelected'] 
#        nQCDantiselVar = bins[njb][stb][htb][btb]['nAntiSelectedVar'] 
##          print nQCDsel, nQCDantisel
#        if nQCDantisel>0:
#          F=nQCDsel/nQCDantisel
#          print 'F_sel-to-anti-sel('+str(stb)+','+str(njb)+','+str(htb)+'):',F
#          if F>0:
#            F_err= F*sqrt(nQCDselVar/nQCDsel**2+nQCDantiselVar/nQCDantisel**2)
#            ratio_ht[stb][njb][btb].SetBinContent(i_htb+1,F)
#            ratio_ht[stb][njb][btb].SetBinError(i_htb+1,F_err)
#            print 'F_sel-to-anti-sel Error('+str(stb)+','+str(njb)+','+str(htb)+'):',F_err
#            bins[njb][stb][htb][btb].update({'F_seltoantiselMC':F, 'F_err':F_err})
#        ratio_ht[stb][njb][btb].GetXaxis().SetBinLabel(i_htb+1, varBinName(htb,'H_{T}'))
#        ratio_ht[stb][njb][btb].GetYaxis().SetTitle('F_{sel-to-antisel}')
#        ratio_ht[stb][njb][btb].GetYaxis().SetRangeUser(0.0,1.0)
##        ratio_ht[stb][njb].GetXaxis().SetTitle('F_{sel-to-antisel}')
#      l.AddEntry(ratio_ht[stb][njb][btb], nJetBinName(njb))
#      if first:
#        ratio_ht[stb][njb][btb].Draw()
#        first = False
#      else:
#        ratio_ht[stb][njb][btb].Draw('same') 
#      l.Draw()
#      t.DrawLatex(0.175,0.85,varBinName(stb,'S_{T}'))
#      text.DrawLatex(0.15,.96,"CMS #bf{#it{Preliminary}}")
#      text.DrawLatex(0.6,0.96,"#bf{L="+str(targetLumi)+" pb^{-1} (13 TeV)}")
#      canv.Print(wwwDir+presel+'Fsa_ht_'+nameAndCut(stb, None, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.png')
#      canv.Print(wwwDir+presel+'Fsa_ht_'+nameAndCut(stb, None, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.pdf')
#      canv.Print(wwwDir+presel+'Fsa_ht_'+nameAndCut(stb, None, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.root')
#
#
##plot F_sel-to-antisel binned in ST for all Njets
#ratio_st={}
#for htb in htreg:
#  ratio_st[htb]={}
#  first = True
#  canv2= ROOT.TCanvas('canv2','canv2',600,600)
#  #canv.SetLogy()
#  l2 = ROOT.TLegend(0.65,0.85,0.95,0.95)
#  l2.SetFillColor(0)
#  l2.SetBorderSize(1)
#  l2.SetShadowColor(ROOT.kWhite)
#  
#  t=ROOT.TLatex()
#  t.SetNDC()
#  t.SetTextSize(0.04)
#  t.SetTextAlign(11)
#  for i_njb, njb in enumerate(njreg):
#    ratio_st[htb][njb]={}
#    for btb in btreg:
#      ratio_st[htb][njb][btb]=ROOT.TH1F('ratio_stHist','ratio_stHist',len(streg),0,len(streg))
#      ratio_st[htb][njb][btb].SetLineColor(ROOT_colors[i_njb])
#      ratio_st[htb][njb][btb].SetLineWidth(2)
#      for i_stb, stb in enumerate(streg):
#        nQCDsel = bins[njb][stb][htb][btb]['nQCDSelected'] 
#        nQCDselVar = bins[njb][stb][htb][btb]['nQCDSelectedVar'] 
#        nQCDantisel = bins[njb][stb][htb][btb]['nAntiSelected'] 
#        nQCDantiselVar = bins[njb][stb][htb][btb]['nAntiSelectedVar'] 
##          print nQCDsel, nQCDantisel
#        if nQCDantisel>0:
#          F=nQCDsel/nQCDantisel
#          print 'F_sel-to-anti-sel('+str(stb)+','+str(njb)+','+str(htb)+'):',F
#          if F>0:
#            F_err= F*sqrt(nQCDselVar/nQCDsel**2+nQCDantiselVar/nQCDantisel**2)
#            ratio_st[htb][njb][btb].SetBinContent(i_stb+1,F)
#            ratio_st[htb][njb][btb].SetBinError(i_stb+1,F_err)
#            print 'F_sel-to-anti-sel Error('+str(stb)+','+str(njb)+','+str(htb)+'):',F_err
#        ratio_st[htb][njb][btb].GetXaxis().SetBinLabel(i_stb+1, varBinName(stb,'S_{T}'))
#        ratio_st[htb][njb][btb].GetYaxis().SetTitle('F_{sel-to-antisel}')
#        ratio_st[htb][njb][btb].GetYaxis().SetRangeUser(0.0,1.0)
##        ratio_st[htb][njb].GetXaxis().SetTitle('F_{sel-to-antisel}')
#      l2.AddEntry(ratio_st[htb][njb][btb], nJetBinName(njb))
#      if first:
#        ratio_st[htb][njb][btb].Draw()
#        first = False
#      else:
#        ratio_st[htb][njb][btb].Draw('same') 
#      l2.Draw()
#      t.DrawLatex(0.175,0.85,varBinName(htb,'H_{T}'))
#      text.DrawLatex(0.15,.96,"CMS Simulation")
#      text.DrawLatex(0.65,0.96,"L="+str(targetLumi/1000)+" fb^{-1} (13 TeV)")
#      canv2.Print(wwwDir+presel+'Fsa_st_'+nameAndCut(None, htb, njetb=None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.png')
#      canv2.Print(wwwDir+presel+'Fsa_st_'+nameAndCut(None, htb, njetb=None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.pdf')
#      canv2.Print(wwwDir+presel+'Fsa_st_'+nameAndCut(None, htb, njetb=None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.root')
#
##plot F_sel-to-antisel binned in ST vs HT
#ratio_2d={}
#for njb in njreg:
#  ratio_2d[njb]={}
#  canv3= ROOT.TCanvas('canv3','canv3',600,600)
#  #canv.SetLogy()
##  l3 = ROOT.TLegend(0.65,0.75,0.95,0.95)
##  l3.SetFillColor(0)
##  l3.SetBorderSize(1)
##  l3.SetShadowColor(ROOT.kWhite)
#  
#  t=ROOT.TLatex()
#  t.SetNDC()
#  t.SetTextSize(0.04)
#  t.SetTextAlign(11)
#  for btb in btreg:
#    ratio_2d[njb][btb]={}
#    ratio_2d[njb][btb]=ROOT.TH2F('ratio_2dHist','ratio_2dHist',len(htreg),0,len(htreg),len(streg),0,len(streg))
#    for i_htb, htb in enumerate(htreg):
#      ratio_2d[njb][btb].GetXaxis().SetBinLabel(i_htb+1,varBinName(htb,'H_{T}'))
#    for i_stb, stb in enumerate(streg):
#      ratio_2d[njb][btb].GetYaxis().SetBinLabel(i_stb+1,varBinName(stb,'S_{T}'))
#
#    for i_htb, htb in enumerate(htreg):
#      for i_stb, stb in enumerate(streg):
#        nQCDsel = bins[njb][stb][htb][btb]['nQCDSelected'] 
#        nQCDselVar = bins[njb][stb][htb][btb]['nQCDSelectedVar'] 
#        nQCDantisel = bins[njb][stb][htb][btb]['nAntiSelected'] 
#        nQCDantiselVar = bins[njb][stb][htb][btb]['nAntiSelectedVar'] 
##          print nQCDsel, nQCDantisel
#        if nQCDantisel>0:
#          F=nQCDsel/nQCDantisel
#          print 'F_sel-to-anti-sel('+str(stb)+','+str(njb)+','+str(htb)+'):',F
#          if F>0:
#            F_err= F*sqrt(nQCDselVar/nQCDsel**2+nQCDantiselVar/nQCDantisel**2)
#            ratio_2d[njb][btb].SetBinContent(i_htb+1,i_stb+1,F)
#            ratio_2d[njb][btb].SetBinError(i_htb+1,i_stb+1,F_err)
#            print 'F_sel-to-anti-sel Error('+str(stb)+','+str(njb)+','+str(htb)+'):',F_err
##      l.AddEntry(ratio_2d[htb][njb][btb], nJetBinName(njb))
#        ratio_2d[njb][btb].Draw('COLZ TEXTE')
#      t.DrawLatex(0.175,0.85,nJetBinName(njb))
#      text.DrawLatex(0.15,.96,"CMS Simulation")
#      text.DrawLatex(0.65,0.96,"L="+str(targetLumi/1000)+" fb^{-1} (13 TeV)") 
#      canv3.Print(wwwDir+presel+'st_vs_ht_'+nameAndCut(None, None, njetb=njb, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.png')
#      canv3.Print(wwwDir+presel+'st_vs_ht_'+nameAndCut(None, None, njetb=njb, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.pdf')
#      canv3.Print(wwwDir+presel+'st_vs_ht_'+nameAndCut(None, None, njetb=njb, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.root')
#
##plot F_sel-to-antisel binned in nJets for all ST bins
#ratio_nj={}
#for htb in htreg:
#  ratio_nj[htb]={}
#  first = True
#  canv4 = ROOT.TCanvas('canv','canv',600,600)
#  #canv.SetLogy()
#  l3 = ROOT.TLegend(0.65,0.80,0.95,0.95)
#  l3.SetFillColor(0)
#  l3.SetBorderSize(1)
#  l3.SetShadowColor(ROOT.kWhite)
#  text = ROOT.TLatex()
#  text.SetNDC()
#  text.SetTextSize(0.04)
#  text.SetTextAlign(11)
#  t3=ROOT.TLatex()
#  t3.SetNDC()
#  t3.SetTextSize(0.04)
#  t3.SetTextAlign(11)
#  for i_stb, stb in enumerate(streg):
#    ratio_nj[htb][stb]={}
#    for btb in btreg:
#      ratio_nj[htb][stb][btb]=ROOT.TH1F('ratio_njHist','ratio_njHist',len(njreg),0,len(njreg))
#      ratio_nj[htb][stb][btb].SetLineColor(ROOT_colors[i_stb])
#      ratio_nj[htb][stb][btb].SetLineWidth(2)
#      for i_njb, njb in enumerate(njreg):
#        nQCDsel = bins[njb][stb][htb][btb]['nQCDSelected'] 
#        nQCDselVar = bins[njb][stb][htb][btb]['nQCDSelectedVar'] 
#        nQCDantisel = bins[njb][stb][htb][btb]['nAntiSelected'] 
#        nQCDantiselVar = bins[njb][stb][htb][btb]['nAntiSelectedVar'] 
##          print nQCDsel, nQCDantisel
#        if nQCDantisel>0:
#          F=nQCDsel/nQCDantisel
#          print 'F_sel-to-anti-sel('+str(stb)+','+str(njb)+','+str(htb)+'):',F
#          if F>0:
#            F_err= F*sqrt(nQCDselVar/nQCDsel**2+nQCDantiselVar/nQCDantisel**2)
#            ratio_nj[htb][stb][btb].SetBinContent(i_njb+1,F)
#            ratio_nj[htb][stb][btb].SetBinError(i_njb+1,F_err)
#            print 'F_sel-to-anti-sel Error('+str(stb)+','+str(njb)+','+str(htb)+'):',F_err
#        ratio_nj[htb][stb][btb].GetXaxis().SetBinLabel(i_njb+1, nJetBinName(njb))
#        ratio_nj[htb][stb][btb].GetYaxis().SetTitle('F_{sel-to-antisel}')
#        ratio_nj[htb][stb][btb].GetYaxis().SetRangeUser(0.0,1.0)
##        ratio_ht[stb][njb].GetXaxis().SetTitle('F_{sel-to-antisel}')
#      l3.AddEntry(ratio_nj[htb][stb][btb], varBinName(stb,'S_{T}'))
#      if first:
#        ratio_nj[htb][stb][btb].Draw()
#        first = False
#      else:
#        ratio_nj[htb][stb][btb].Draw('same') 
#      l3.Draw()
#      t3.DrawLatex(0.2,0.85,varBinName(htb,'H_{T}'))
#      text.DrawLatex(0.15,.96,"CMS Simulation")
#      text.DrawLatex(0.65,0.96,"L="+str(targetLumi/1000)+" fb^{-1} (13 TeV)")
#      canv4.Print(wwwDir+presel+'Fsa_nj_'+nameAndCut(None, htb, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.png')
#      canv4.Print(wwwDir+presel+'Fsa_nj_'+nameAndCut(None, htb, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.pdf')
#      canv4.Print(wwwDir+presel+'Fsa_nj_'+nameAndCut(None, htb, None, btb=btb, presel="(1)", charge="", btagVar = 'nBJetMediumCSV30')[0]+'.root')


