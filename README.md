# SiPMCalib Analysis Software

The analysis software for the data collected in the [SiPMCalib control
software][sipmcalibcontrol] to extract both SiPM parameters and control
performance metrics. Summary of results of the analysis will be displayed on the
[offical documentation twiki page][sipmcalibtwiki].

For installation, you will need a CMSSW environment that is compatible with C++17
standards, currently tested with `CMSSW_10_3_1_patch1`. Anything newer should
also work. The CMSSW environment is only used for easier code organization.

## Installation on CMSSW comaptible machine

On a machine with that support the correct architecture, run the following commands:

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_3_1_patch1
cd CMSSW_10_3_1_patch1/src
cmsenv

git clone --recurse-submodules https://github.com/yimuchen/UserUtils.git
git clone https://github.com/yimuchen/SiPMCalib.git

scram b
```

This will give you the plotting commands in
`CMSSW_10_3_1_patch1/bin/slc_amd64_gcc700/`

For the documentation of each of the plotting commands, run the plotting command
with the help flag: `<cmd> --help`, which will generate list of the command options.

For the detailed documentation of each sub-module. Look at the `README.md` file in each
of the subdirectories:

- [Common](Common): Common parsing the and formatting routines.
- [InvSqCalc](InvSqCalc): The standard plotting routine for inverse square-law
  related calculation routines.
- [Pulser](Pulser): LED Pulser board parameter extractions
- [SiPMCalc](SiPMCalc): SiPM parameter extraction routines.
- [VisCalc](VisCalc): Standard plotting routines for the visual calibration
  routines.
- [TimeCalc](TimeCalc): Standard analysis routines for temporal scans for
  stability calibration.

Standardized analysis and plotting routines to be called by the user should be
kept in the `bin` and `src` directories of their respective subdirectories.
Testing scripts for quick plots should be kept in the testing directory.

[sipmcalibcontrol]: https://github.com/yimuchen/SiPMCalibControl
[sipmcalibtwiki]: https://twiki.cern.ch/twiki/bin/viewauth/CMS/UMDHGCalSiPMCalib
