# SiPMCalib Analysis Software

The analysis software for the data collected in the [SiPMCalib control
software][SiPMCalibControl] to extract both SiPM parameters and control
performance metrics. Summary of results of the analysis will be displayed on the
[offical documentation twiki page][SiPMCalibTwiki].

For installation, you will need a CMSSW environment that is compatible with C++17
standards, currently tested with `CMSSW_10_3_1_patch1`. Anything newer should
also work. The CMSSW environment is only used for easier code organization. A
standalone version should be coming soon.

## Installation on CMSSW comaptible machine

On a machine with that support the correct architecture, run the following commands:

```bash
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_3_1_patch1
cd CMSSW_10_3_1_patch1/src
cmsenv

git clone https://github.com/yimuchen/UserUtils.git
git clone https://github.com/yimuchen/SiPMCalib.git

scram b
```

This will give you the plotting commands in
`CMSSW_10_3_1_patch1/bin/slc_amd64_gcc700/`

For the documentation of each of the plotting commands, check the individual
documentation in each of the subdirectories. Or run the `<cmd> --help` to pull up
the list of the command options.

[SiPMCalibControl]: https://github.com/yimuchen/SiPMCalibControl
[SiPMCalibTwiki]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/UMDHGCalSiPMCalib
