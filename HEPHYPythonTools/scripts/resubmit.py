#!/bin/env python
import os, sys, re
infile=""
if len(sys.argv)==1:
  print "No input file"
#  infile = "status.txt"
  sys.exit()
else:
  infile = sys.argv[1]
print "Input File", infile
lines = file(infile).readlines()
resubmitjobs = {}
latestwdir = "None"
wdirPat = re.compile('	working directory')
joblinePat = re.compile('^[0-9][0-9]* ')

indications = [\
"Aborted",
"7002",
"8020",
"8004",
"8009",
"8021",
"8018",
"8028",
"8001",
"8002",
"10016",
"10034",
"50117",
"60317",
"60307",
"60302",
"60318",
"50115",
"50669",
"50664",
"70500",
"50700",
"50800",
"Cancelled",
"1          1",
"9          9",
"127        127",
"Created"
]
for line in lines:
  if wdirPat.match(line):# or lines.count
    wstring = line.replace("working directory", "").replace('	','').replace(" ","").replace('\n','').replace('\t','')
    resubmitjobs[wstring] = []                        
    latestwdir = wstring                              
  jobm = joblinePat.match(line)                       
  if jobm:                                            
    for ind in indications:                           
      if line.count(ind)>0:                           
        jobn  = int(jobm.group())
        resubmitjobs[latestwdir] .append(jobn)
        break
outfile = file('resubmit.sh', 'w')
outfile.write('#!/bin/sh\n')
for key in resubmitjobs:
  if len(resubmitjobs[key])>0:
    prestring = "crab -forceResubmit "
    for n in range(len(resubmitjobs[key])/500+1):
      sstring = prestring + ','.join([str(j) for j in resubmitjobs[key][500*n:500*(n+1)]])+' -c ' +key
      outfile.write(sstring+'\n') 
      print sstring
print "Written resubmit.sh"
outfile.close()
os.system("chmod +x resubmit.sh")