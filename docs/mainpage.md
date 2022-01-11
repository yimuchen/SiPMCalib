# SiPM calibration analysis package

Analysis package for the SiPM calibration project at the UMD CMS group. This
package is intended to work in conjunction with the project's data collection
program (found [here][sipmcalibcontrol]), analyzing the data and parsing the data
into a publication ready plot.

## Installation

This package should be deployed in a CMSSW environment, but doesn't actually
depend on any CMSSW packages. The choice of using CMSSW is purely to help with
the simplification of the code organization and compiling process for C++ code.
This package also relies on the package [UserUtils][userutils] for the generation
of publication ready plots.

Any `CMSSW` version/`SCRAM_ARCH` compatible with C++17 would allow the package to
be compiled, though the tested most thoroughly with `CMSSW_10_3_X` (purely
historical reason).

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




[sipmcalibcontrol]: https://github.com/umdcms/SiPMCalibControl
[userutils]: https://github.com/yimuchen/UserUtils
