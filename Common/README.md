# SiPM calibration software --- Common

Here contains the common functions used for parsing the standardized readout
formats that is generated using the data collection routines in the SiPM
calibration control software.

---

## [StdFormat](interface/StdFormat.hpp)

The standard format for general data collection used by the gantry system. The
data will contain the following space separated columns in an ASCII file:

- The time stamp (`double`, seconds)
- The detector ID (`int`)
- The gantry coordinate of data collection (`double, double, double`, mm)
- The bias voltage of the LED pulser (`double`, mV)
- The temperature of the LED pulser (`double`, C)
- The temperature of the SiPM board (`double`, C)

All columns after the 8th column is regarded as `data`, the meaning of which will
depend on the data collection process of interest. Is the most common case, the
data itself will only contain `2N` columns representing some data readout and
uncertainty of the readout respectively. Here, we would typically want to extract
each column into a double array of the same length for additional data processing.
In this case, one would typically use the methods like:

```c++
StdFormat::X()
StdFormat::DataCol(1)
```

To extract columns of data into a `std::vector` container. A function pointer can
be passed as the second argument to select only particular rows of interest.

---

## [WaveFormat](interface/SaveFormat.hpp)

The standard format for waveform like data collection used by the gantry system.
The first row needs to contain the following 3 numbers:

- The sampling time for each sample (`double`, ns).
- The number of bits used for a sample (`int`).
- The ADC convertion rate for a single bit (`double`, mV).

The remaining rows should container `2xNxM` alpha-numerical numbers, representing
M sample measurements of a single waveform in N-bit hexadecimal format. Once the
class processes the input file, the user can then use the inbuilt function to run
routines like summing the waveform with or without pedestal subtraction and such.
