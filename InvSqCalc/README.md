# SiPMCalib Analysis Software --- InvSqCalc

Plotting and calculations scripts regarding the performance on inverse-square law
performance. This package makes no assumption on the type of light-source
(dynamic or static), as long as the file is in the expected format.

## InvSq_UpdateMCFile

Running a Monte-Carlo calculation for the expected light yield received at a
rectangular sensor by a small, disk-like light source. The user can specify a
list of aperture diameters, vertical distance, horizontal displacements and SiPM
size. The format will be stored as a plain text file, with the pass/fail of
randomly generated rays stored in the file.

## InvSq_PlotMC

Plot a list of MC results in a Luminosity--z position format.

## InvSq_HScanPlot

Plot the results of a luminosity scan in the x-y plane of the gantry coordinates,
along with a model fit prediction using an ideal inverse square law. The stuff
should be similar to the numpy results obtained during the gantry run.

## InvSq_ZScanPlot

Plot the results of a luminosity scan in the z direction of the gantry coordinate
system, along with a model fit prediction using an ideal inverse square law with
potentially non-zero z-offset, horizontal offset and pedestal. The plot will be
shifted such that z-offset, horizontal offset is consistent with zero, for easier interpretation of the fit result.

## InvSq_CompZProfile

Given two luminosity scan data files, plot the luminosity--z profile for both files, the second luminosity will be normalized to the first.
