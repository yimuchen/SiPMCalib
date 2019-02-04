#!/bin/env python3

template = """
universe = vanilla
Executable = /data/users/yichen/CMSSW_10_3_1_patch1/src/SiPMCalib/SiPMCalc/condor_runtoy.sh
should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stdout
Error  = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stderr
Log    = /data/users/yichen/condor/log/simple_$(cluster)_$(process).condor
Arguments = --nToys 1000 --nEvents {0} --model {1} --fit {2} --rate {3} --xtalkrate {4} --darkrate {5} alpha {6}
Queue 1
"""

nEvents = [3000,10000,30000,100000]
models  = ["full","ndc","nac","simp","dark"]
fits    = ["binned","unbinned"]
means   = [1.5,4,12]
xrates  = [0.1,0.2,0.3]
dcfracs = [0.05,0.1,0.2]
alphas  = [0.05,0.1,0.2]

for nEvent,model,fit,mean,xrate,dcfrac,alpha in [
    (a,b,c,d,e,f,g) for a in nEvents 
                      for b in models
                      for c in fits
                      for d in means
                      for e in xrates
                      for f in dcfracs 
                      for g in alphas  ]:
  filename = "/data/users/yichen/condor/jdl/runtoy_{0}_{1}_{2}_{3}_{4}_{5}_{6}.jdl".format(nEvent,model,fit,mean,xrate,dcfrac,alpha)
  
  print(filename)
  file = open(filename, "w")
  file.write(template.format(nEvents,model,fit,mean,xrate,dcfrac,alpha))
