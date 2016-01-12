#!/bin/sh 
########Spring15###############

#python cmgPostProcessing.py --leptonSelection=hard --skim=""  --samples=TToLeptons_sch_25ns
#python cmgPostProcessing.py --leptonSelection=hard --skim=""  --samples=TToLeptons_tch_25ns
#python cmgPostProcessing.py --leptonSelection=hard --skim=""  --samples=TBar_tWch_25ns
#python cmgPostProcessing.py --leptonSelection=hard --skim=""  --samples=T_tWch_25ns

#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TToLeptons_sch
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TToLeptons_tch
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TBar_tWch
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=T_tWch
#
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TTWJetsToLNu_25ns
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TTWJetsToQQ_25ns
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TTZToLLNuNu_25ns
#python cmgPostProcessingAntiSelectionV2.py --overwrite --leptonSelection=none --samples=TTZToQQ_25ns

python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TToLeptons_sch
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TToLeptons_tch
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TBar_tWch
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=T_tWch

python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TTWJetsToLNu_25ns
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TTWJetsToQQ_25ns
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TTZToLLNuNu_25ns
python cmgPostProcessing.py --overwrite  --skim="HT500ST250"  --samples=TTZToQQ_25ns

