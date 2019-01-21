#!/bin/env python
template = """
universe = vanilla
Executable = /data/users/yichen/CMSSW_10_3_1_patch1/src/SiPMCalib/SiPMCalc/condor_runtoy.sh
should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stdout
Error  = /data/users/yichen/condor/log/simple_$(cluster)_$(process).stderr
Log    = /data/users/yichen/condor/log/simple_$(cluster)_$(process).condor
Arguments = --nToys 1000 --nEvents {0} --model {1}
Queue 1
"""

for nEvents in [3000,10000,30000,100000] :
  for model in ["full","gpgaus", "pgaus"]:
    filename = "/data/users/yichen/condor/jdl/runtoy_{0}_{1}.jdl".format(model,nEvents)
    print(filename)
    file = open(filename, "w")
    file.write(template.format(nEvents,model))
