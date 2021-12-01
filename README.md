# SiPMCalib Analysis Software

The analysis software for the data collected in the [SiPMCalib control
software][sipmcalibcontrol] to extract both SiPM parameters and control
performance metrics. Summary of results of the analysis will be displayed on the
[offical documentation twiki page][sipmcalibtwiki]. The main documentation for
this package can be found [here][maindoc].

## Installation on CMSSW compatible machine

On a machine with that support the correct architecture, run the following commands:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_3_1_patch1
cd CMSSW_10_3_1_patch1/src
cmsenv

git clone --recurse-submodules https://github.com/yimuchen/UserUtils.git
git clone https://github.com/yimuchen/SiPMCalib.git

scram b
```

[sipmcalibcontrol]: https://github.com/yimuchen/SiPMCalibControl
[sipmcalibtwiki]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/UMDHGCalSiPMCalib
[maindoc]: https://yimuchen.github.io/SiPMCalib/
