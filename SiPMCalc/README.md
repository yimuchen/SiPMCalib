# SiPMCalib Analysis Software -- SiPMCalc

The functions provided in this function provides for the extraction of the
various SiPM operation parameters. Sources of the analysis methods will be given
in the documentation of each of the binaries.

---

## SiPM_DarkTrigger

Given the waveform file of the SiPM readout triggered via the random dark current
discharges, compute the secondary noise parameters: cross-talk and after pulsing
probabilty and time-scale, as well as the SiPM recovery time.

The program takes just two arguments, the input file and the prefix of the output
file. Computation results will be printed on screen, as well as printed in the
generated plots.

---

## SiPM_FitDark

Given the waveform file of a SiPM readout that is triggered via a random
trigger in a no-light condition, compute the primary and secondary noise
parameters of the SiPM using the area spectrum of waveforms.

In the program, one can specify where the integration window should start and
end, as well as the binning scheme of the data. It is advised to relatively wide
binning in the data, the data is surprisingly discrete in both the voltage and
temporal resolution.

---

## SiPM_FitLowLight

Given either a waveform file or a the standard data file using the
lowlightcollect command of a SiPM readout collected by the gantry system, compute
all operation parameters characterizing a SiPM low light response. The user will
be reponsible for ensuring that the save files is provided.

In the program, one can specify where the integration window should start and
end, the binning scheme of data to be used. All parameters of the PDF functions
can be fixed or given a custom range. Estimations for the pedestal, gain,
Gaussian noise and number of photon will be estimated using a peak finding
algorithm is an estimation was not given.

---

## SiPM_DisplayWaveform

Given a waveform file of a SiPM readout, display the waveform traces as a heat map.
