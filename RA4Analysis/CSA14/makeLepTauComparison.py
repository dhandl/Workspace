import ROOT
from Workspace.RA4Analysis.stage2Tuples import ttJetsCSA14

htCut     = 400
metCut    = 150
minNJets  =   4

#htCut     = 0
#metCut    = 0
#minNJets  =   0

c = ROOT.TChain('Events')
for b in ttJetsCSA14['bins']:
  c.Add(ttJetsCSA14['dirname']+'/'+b+'/h*.root')

cPred = ROOT.TChain('Events')
cPred.Add('/data/schoef/results2014/tauTuples/CSA14_TTJets_lepGenTau.root')

#oneLepTau     ="ngoodMuons==1&&nvetoMuons==1&&nvetoElectrons==0&&Sum$(gTauPt>15&&abs(gTauEta)<2.5&&gTauNENu+gTauNMuNu==1&&gTauNTauNu==1)==1"
#oneLepTauOpen ="ngoodMuons==1&&nvetoMuons==1&&nvetoElectrons==0&&Sum$(gTauNENu+gTauNMuNu==1&&gTauNTauNu==1)==1"
oneLepTau     ="ngNuMuFromW==1&&ngNuEFromW==0&&ngNuTauFromW==1&&ngoodMuons==1&&nvetoMuons==1&&nvetoElectrons==0&&Sum$(gTauPt>15&&abs(gTauEta)<2.1&&gTauNENu==0&&gTauNMuNu==1&&gTauNTauNu==1)==1"
oneLepTauOpen ="ngNuMuFromW==1&&ngNuEFromW==0&&ngNuTauFromW==1&&ngoodMuons==1&&nvetoMuons==1&&nvetoElectrons==0&&Sum$(gTauNENu==0&&gTauNMuNu==1&&gTauNTauNu==1)==1"

#scaleF = (1-0.1741-0.1783)*0.1125/(0.1057+0.1075)
#scaleF = (0.1741+0.1783)*0.1125/(0.1057)
scaleF = (0.1741)*0.1125/(0.1057)

hMT = ROOT.TH1F('hMT', 'hMT',20,0,800)
c.Draw('sqrt(2.*met*leptonPt*(1-cos(leptonPhi-metphi)))>>hMT','weight*('+oneLepTau+'&&ht>'+str(htCut)+'&&njets>='+str(minNJets)+'&&met>'+str(metCut)+')','goff')

hMTOpen = ROOT.TH1F('hMTOpen', 'hMTOpen',20,0,800)
c.Draw('sqrt(2.*met*leptonPt*(1-cos(leptonPhi-metphi)))>>hMTOpen','weight*('+oneLepTauOpen+'&&ht>'+str(htCut)+'&&njets>='+str(minNJets)+'&&met>'+str(metCut)+')','goff')

hMTPred = ROOT.TH1F('hMTPred', 'hMTPred',20,0,800)
cPred.Draw('mTPred>>hMTPred','weightPred*scaleLEff*(htPred>'+str(htCut)+'&&njetsPred>='+str(minNJets)+'&&metPred>'+str(metCut)+')','goff')

c1=ROOT.TCanvas()
hMT.Draw()
hMTPred.SetLineColor(ROOT.kRed)
hMTPred.Scale(scaleF)
#hMTPred.Scale(hMT.Integral()/hMTPred.Integral())
hMTPred.Draw('same')

hMTOpen.SetLineStyle(2)
hMTOpen.Draw('same')

c1.SetLogy()
c1.Print('/afs/hephy.at/user/s/schoefbeck/www/pngCSA14/compMT_lepTau.png')

hHT = ROOT.TH1F('hHT', 'hHT',25,0,2500)
c.Draw('ht>>hHT','weight*('+oneLepTau+'&&ht>'+str(htCut)+'&&njets>='+str(minNJets)+'&&met>'+str(metCut)+')','goff')

hHTPred = ROOT.TH1F('hHTPred', 'hHTPred',25,0,2500)
cPred.Draw('htPred>>hHTPred','weightPred*scaleLEff*(htPred>'+str(htCut)+'&&njetsPred>='+str(minNJets)+'&&metPred>'+str(metCut)+')','goff')

c1=ROOT.TCanvas()
hHTPred.SetLineColor(ROOT.kRed)
hHTPred.Scale(scaleF)
#hHTPred.Scale(hHT.Integral()/hHTPred.Integral())
hHTPred.Draw()
hHT.Draw('same')
c1.SetLogy()
c1.Print('/afs/hephy.at/user/s/schoefbeck/www/pngCSA14/compHT_lepTau.png')

hNJet = ROOT.TH1F('hNJet', 'hNJet',12,0,12)
c.Draw('njets>>hNJet','weight*('+oneLepTau+'&&ht>'+str(htCut)+'&&njets>='+str(minNJets)+'&&met>'+str(metCut)+')','goff')

hNJetPred = ROOT.TH1F('hNJetPred', 'hNJetPred',12,0,12)
cPred.Draw('njetsPred>>hNJetPred','weightPred*scaleLEff*(htPred>'+str(htCut)+'&&njetsPred>='+str(minNJets)+'&&metPred>'+str(metCut)+')','goff')

c1=ROOT.TCanvas()
hNJetPred.SetLineColor(ROOT.kRed)
hNJetPred.Scale(scaleF)
#hNJetPred.Scale(hNJet.Integral()/hNJetPred.Integral())
hNJetPred.Draw()
hNJet.Draw('same')
c1.SetLogy()
c1.Print('/afs/hephy.at/user/s/schoefbeck/www/pngCSA14/compNJet_lepTau.png')
