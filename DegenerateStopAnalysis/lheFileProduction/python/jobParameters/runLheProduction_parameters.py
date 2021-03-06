#
# default parameters for LHE production
# 
# parameters can be overwritten for a given job 

import os

class LheProductionParameters:


    def __init__(self, jobParametersFile):

        self.stopMassLimits = getattr(__import__(jobParametersFile, 
                           fromlist=['stopMassLimits']), 'stopMassLimits')
        
        self.generatedLspMassLimits = getattr(__import__(jobParametersFile, 
                           fromlist=['generatedLspMassLimits']), 'generatedLspMassLimits')

        self.deltaMassStopLsp = getattr(__import__(jobParametersFile, 
                           fromlist=['deltaMassStopLsp']), 'deltaMassStopLsp')

        self.splitUndecayedSample = getattr(__import__(jobParametersFile, 
                           fromlist=['splitUndecayedSample']), 'splitUndecayedSample')

        self.deltaMassStopLspFractions = getattr(__import__(jobParametersFile, 
                           fromlist=['deltaMassStopLspFractions']), 'deltaMassStopLspFractions')

        self.deltaMassStopLspSelected = getattr(__import__(jobParametersFile, 
                           fromlist=['deltaMassStopLspSelected']), 'deltaMassStopLspSelected')
        
        if (len(self.deltaMassStopLspSelected) == 0):
            self.deltaMassStopLspSelected = self.deltaMassStopLsp

        self.workDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['workDirectory']), 'workDirectory')

        self.madGraphDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['madGraphDirectory']), 'madGraphDirectory')

        self.undecayedFilesDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['undecayedFilesDirectory']), 'undecayedFilesDirectory')

        self.mergedFilesDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['mergedFilesDirectory']), 'mergedFilesDirectory')

        self.lheUndecayedEosDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['lheUndecayedEosDirectory']), 'lheUndecayedEosDirectory')

        self.undecayedFilesStageDirectory = getattr(__import__(jobParametersFile, 
                           fromlist=['undecayedFilesStageDirectory']), 'undecayedFilesStageDirectory')

        self.stageUndecayedLheFiles = getattr(__import__(jobParametersFile, 
                           fromlist=['stageUndecayedLheFiles']), 'stageUndecayedLheFiles')

        self.numberEvents = getattr(__import__(jobParametersFile, 
                           fromlist=['numberEvents']), 'numberEvents')
    
    def __str__(self):

        lheUndecayedEosDirectoryStr =     '    EOS LHE undecayed sample:       ' + self.lheUndecayedEosDirectory
        undecayedFilesStageDirectoryStr = '    Staged LHE undecayed sample:    ' + self.undecayedFilesStageDirectory
        undecayedFilesDirectoryStr =      '    Job local LHE undecayed sample: ' + self.undecayedFilesDirectory
        mergedFilesDirectoryStr =         '    Merged LHE decayed sample:      ' + self.mergedFilesDirectory
        workDirectoryStr =                '    Work directory:                 ' + self.workDirectory
        
        if self.stopMassLimits[0] < 0:
            stopMassLimitsStr = '    Use all available stop mass values'
        else:
            stopMassLimitsStr = '    ' + str(self.stopMassLimits[0]) + ' <= stop mass <= ' + str(self.stopMassLimits[1])
    
        if self.generatedLspMassLimits[0] < 0:
            generatedLspMassLimitsStr = '    Use all generated LSP mass values'
        else:
            generatedLspMassLimitsStr = '    ' + str(self.generatedLspMassLimits[0]) + ' <= generated LSP mass <= ' + str(self.generatedLspMassLimits[1])
                         
        deltaMassStopLspStr = ''
        for dMassIndex, dMass in enumerate(self.deltaMassStopLspSelected):
            if (self.splitUndecayedSample == True):
                deltaMassStopLspStr += '        deltaMassStopLsp = ' + str(dMass) + ' sample fraction =' +  str(self.deltaMassStopLspFractions[dMassIndex]) + '\n'
            else:
                deltaMassStopLspStr += '        deltaMassStopLsp = ' + str(dMass) + '\n'
     
        if (self.splitUndecayedSample == True):
            splitUndecayedSampleStr = '    splitUndecayedSample = True'
        else:
            splitUndecayedSampleStr = '    splitUndecayedSample = False'
     
        if (self.stageUndecayedLheFiles == True):
            stageUndecayedLheFilesStr = '    stageUndecayedLheFiles = True'
        else:
            stageUndecayedLheFilesStr = '    stageUndecayedLheFiles = False'

        if self.numberEvents < 0:
            numberEventsStr = '    Use all available events'
        else:
            numberEventsStr = '    Number of decayed events = ' + str(self.numberEvents)

        return \
            '\nInput parameters' + '\n' + \
            lheUndecayedEosDirectoryStr + '\n' + \
            undecayedFilesStageDirectoryStr + '\n' + \
            workDirectoryStr + '\n' + \
            undecayedFilesDirectoryStr + '\n' + \
            mergedFilesDirectoryStr + '\n' + \
            stopMassLimitsStr + '\n' + \
            generatedLspMassLimitsStr + '\n' + \
            '    deltaMassStopLsp values' + '\n' + \
            deltaMassStopLspStr + '\n' + \
            splitUndecayedSampleStr + '\n' + \
            stageUndecayedLheFilesStr + '\n' + \
            numberEventsStr + '\n' + \
            '\n'
        
    def checkConsistency(self):
        # check input parameters for consistency

        for dMass in self.deltaMassStopLspSelected:
            if dMass not in self.deltaMassStopLsp:
                sys.exit('\nInconsistent deltaMassStopLspSelected and deltaMassStopLsp parameters. Selected values are subset of full list')

        if ((len(self.deltaMassStopLsp) != len(self.deltaMassStopLspFractions)) and  (self.splitUndecayedSample == True)):
            sys.exit('\nInconsistent deltaMassStopLsp and deltaMassStopLspFractions parameters. Different sizes.')
      
        if ((self.undecayedFilesStageDirectory) == ''):
            sys.exit('\nEmpty directory name for the stage directory of undecayed LHE file.')

        if ((self.undecayedFilesDirectory) == ''):
            sys.exit('\nEmpty directory name for the job local directory of undecayed LHE file.')

        if ((self.mergedFilesDirectory) == ''):
            sys.exit('\nEmpty directory name for the job local directory of merged LHE file.')

# the stop mass values are fixed for the LHE files with undecayed stops, 
# one can only choose the range the files are selected and processed
#
# choose the range the files are selected and processed
# if stopMassLimits[0] is negative, it will process files for all existing stop mass values
# otherwise it will process all files having 
#     stopMassLimits[0] <= stop mass <= stopMassLimits[1]

# generated stop mass limits (GeV)
stopMassLimits = [-1]
#stopMassLimits = [0, 1000]

# there are no LSP particles in the undecayed files (stop and antistop are not decayed)
# the generated LSP mass values are, as such, dummy, but they are included in the file name
#
# choose the generated LSP mass values for which the files are selected and processed by MadGraph
# if generatedLspMass[0] is negative, it will process all the files corresponding to the stop mass(es) given above
# otherwise it will process all files having 
#     generatedLspMassLimits[0] <= generated LSP mass <= generatedLspMassLimits[1]
generatedLspMassLimits = [-1] 
#generatedLspMassLimits = [50, 75]

# mass difference stop - LSP (GeV)
# one can choose the mass difference stop - LSP freely (stop and antistop are not decayed)
# deltaMassStopLsp - list of all mass differences to be considered
# splitUndecayedSample - if True, one undecayed file is used for one deltaMassStopLsp value only
# deltaMassStopLspFractions - fraction in percents of the undecayed events to be used for a given deltaMassStopLsp
#    Note: files are not split, so the fraction is approximate 
#
# deltaMassStopLspSelected - list of mass differences selected from deltaMassStopLsp to be processed in this job
# if deltaMassStopLspSelected has no elements, it uses the values from deltaMassStopLsp

deltaMassStopLsp = [10, 20, 30, 40, 50, 60, 70, 80]
splitUndecayedSample = False
deltaMassStopLspFractions = [0.125]*8

deltaMassStopLspSelected = []

# number of decay events to be generated
# if negative, take the number from  official_run_card_decay_TEMPLATE.dat
numberEvents = -1
#numberEvents = 100

# work directory for a job
workDirectory = os.getcwd()

# actual values for the empty names of the directories will be set in 
# the run shell script 

# MadGraph home
madGraphDirectory = ''

# directory where undecayed LHE file are to be found by a job 
undecayedFilesDirectory = ''

# directory to save the merged files by a job
mergedFilesDirectory = ''

# EOS location of the T2tt undecayed LHE files
# 
# lheSample: T2tt or stop_stop

lheSample = 'stop_stop'

if lheSample == 'stop_stop':
    lheUndecayedEosDirectory = '/store/group/phys_susy/LHE/stop_stop/T2tt_Undecayed'
elif lheSample == 'T2tt':   
    lheUndecayedEosDirectory = '/store/group/phys_susy/LHE/T2tt'
else:
    sys.exit('No valid EOS LHE sample.')

# directory to stage undecayed LHE file, or where files are already staged
# directory must be given as absolute path
undecayedFilesStageDirectory = ''

# options to control the job
#

# stageUndecayedLheFiles
#     if the value is False, the program assumes the files are already staged in a local/storage element directory 
#         defined by undecayedFilesStageDirectory
#         ! in this case, the list of files to be processed is built from that directory (only)
#     if True, it stages the selected files from EOS directory defined by lheUndecayedEosDirectory, then exit
stageUndecayedLheFiles = False



