#!/bin/env python3

template = """
universe = vanilla
Executable = /data/users/yichen/CMSSW_10_3_1_patch1/src/SiPMCalib/SiPMCalc/condor_runtoy.sh
should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stdout
Error  = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stderr
Log    = /data/users/yichen/condor/log/simple_$(cluster)_$(process).condor
Arguments = --nToys 500 --nEvents {0} --model {1} --fit {2} --rate {3} --xtalkrate {4} --darkrate {5} alpha {6}
Queue 1
"""

nEvents = [3000,10000,30000,100000]
models  = ["full","ndc","nap","simp","dark"]
fits    = ["binned","unbinned"]
means   = [1.5,4,12]
xrates  = [0.1,0.2,0.3]
dcfracs = [0.05,0.1,0.2]
alphas  = [0.05,0.1,0.2]

for nEvent,model,fit,mean,xrate,dcfrac,alpha in [
    (a,b,c,d,e,f,g) for a in nEvents
                      for b in models                      for c in fits
                      for d in means
                      for e in xrates
                      for f in dcfracs
                      for g in alphas  ]:

  filename = "SiPMRunToy_{0}_nEvt{1}_{2}".format(model,nEvent,fit)
  if model != 'dark':
    filename = filename + "_r{0}_x{1}".format(mean,xrate)
  if model != 'ndc' and model != 'simp':
    filename = filename + "_dc{0}".format(dcfrac)
  if model != 'nap' and model != 'simp':
    filename =filename + "_a{0}".format(alpha)
  filename = "/data/users/yichen/condor/jdl/" + filename.replace('.','p') + ".jdl"

  print(filename)
  file = open(filename, "w")
  file.write(template.format(nEvent,model,fit,mean,xrate,dcfrac,alpha))
