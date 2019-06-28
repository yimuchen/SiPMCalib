#!/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/yichen/local/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/yichen/local/lib64/
cd /data/users/yichen/CMSSW_10_3_1_patch1/src/SiPMCalib/SiPMCalc
eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
/data/users/yichen/CMSSW_10_3_1_patch1/bin/slc6_amd64_gcc700/SiPM_RunToy $@
