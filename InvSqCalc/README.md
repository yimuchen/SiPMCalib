# SiPMCalib Analysis Software --- InvSqCalc

Plotting and calculations scripts regarding the performance on inverse-square law
performance. This package makes no assumption on the type of light-source
(dynamic or static), as long as the file is in the expected format.

## InvSq_HScanPlot

Given a data file representing the results of a luminosity scan in the x-y plane
of the gantry, fit and plot the results of comparing the data with a simple
inverse-square law. This is rarely used for calibration purposes, and many used
to generate plots used in notes and papers.

## InvSq_ZScanPlot

Given a data file representing the results of a luminosity scan in the z
direction of the gantry, fit and plot the results of comparing the data with a
simple inverse-square law. The plot will be present the data with shifted
z-offsets and pedestal value according to the estimate obtained in fit estimation
for easier interpretation. If the photo-detector is linear (a.k.a a photo diode),
the fitted profile can be used as a reference for the non-linearity fitting.
