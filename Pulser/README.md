# SiPMCalib Analysis Software --- PulseShape

Plotting commands for extracting the performance of the LED pulser board.

## SiPM_PulseShape

Take the input the Texas Instrument TDC7201 evaluation board and use it to plot
the leading edge distribution of the LED light source. User is responsible for
specifying the range of and binning of the plot to get a reasonable result. The
input file format assumes that the leading 4 rows contains TDC7201 header
information (which is discarded), while only the first columns of the remaining
rows will be used in the plot.
