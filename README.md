# Detector Spectrum Analysis (ROOT)

ROOT macros for reading, plotting, and fitting a charged-particle detector energy spectrum. The spectrum contains 5 overlapping peaks sitting on a background; the macros fit them with a sum of Gaussians plus a linear or exponential background, and report the position, area, and uncertainty of each peak.

Two variants are provided, differing only in what the histogram's X-axis represents:

| File | X-axis | Use case |
|---|---|---|
| `AnalyzeDetectorSpectrum_Channel.C` | Raw ADC channel | Fitting is done directly on raw channels; an energy calibration is used only to *label* each fitted peak's energy (keV) in the printed output. |
| `AnalyzeDetectorSpectrum_Calibrated.C` | Calibrated energy (keV) | The histogram itself is built with a calibrated energy axis, so fit inputs (ROI, peak guesses, sigma guesses) are given directly in keV. |

## Repository contents

```
.
├── AnalyzeDetectorSpectrum_Channel.C      # fit on raw ADC channel axis
├── AnalyzeDetectorSpectrum_Calibrated.C   # fit on calibrated energy (keV) axis
├── Histo_test.txt                        # raw spectrum data (input)
└── README.md
```

`Histo_test.txt` is a plain text file with one integer per line: the number of counts in that ADC channel. There is no header; line 1 = channel 0, line 2 = channel 1, and so on.

## Requirements

- [ROOT](https://root.cern) 6.x (tested with the standard `TH1F`, `TF1`, `TGraphErrors`, `TFitResult` classes; no external dependencies beyond a standard ROOT installation).

## Usage

Run either macro directly from the ROOT prompt or the command line. ROOT automatically calls the function that matches the file name.

```bash
# Full analysis: fit background + 5 peaks, print results
root -l AnalyzeDetectorSpectrum_Channel.C
root -l AnalyzeDetectorSpectrum_Calibrated.C

# Display-only mode: just read and plot the spectrum, skip fitting entirely
root -l 'AnalyzeDetectorSpectrum_Channel.C(true)'
root -l 'AnalyzeDetectorSpectrum_Calibrated.C(true)'
```

Both macros expect `Histo_test.txt` to be in the current working directory (or edit the `INPUT_FILE` variable near the top of the `main` function).

### Output

- A canvas (`detector_spectrum.png`) with two pads:
  - **Pad 1** – the full spectrum (log Y), with the region of interest (ROI) marked by two dashed vertical lines.
  - **Pad 2** – a zoomed view of the ROI (log Y), showing the data points, the individual fitted peaks, the fitted background, and the total fit curve.
- Fit results printed to the console: background parameters, per-peak `Amp`/`Mean`/`Sigma` with uncertainties and correlations, overall `Chi2/NDF`, and each peak's integrated area with propagated uncertainty (and its calibrated energy, when available).

## What each macro does

1. **Energy calibration** – fits a linear function `E(keV) = a·channel + b` to 4 known (channel, energy) calibration points.
2. **Read the spectrum** – loads `Histo_test.txt` into a `TH1F` (bin *i+1* = channel *i*; bin edges shifted by `-0.5`/`n-0.5` so bin centers align exactly with the calibrated channel energies).
3. **Background estimation** (`FitBackgroundSides`) – fits a linear (`pol1`) or exponential (`expo`) background using *only* two user-defined peak-free "side windows", one below and one above the peak region, so the peaks themselves cannot bias the background estimate.
4. **Multi-peak fit** (`FitGaussPeaks`) – fits the sum of `N_PEAKS` Gaussians plus the background (background parameters are fixed from step 3, only the peak parameters are free) in a single combined fit.
5. **Peak area calculation** (`ComputePeakAreas`) – computes each Gaussian's integrated area analytically (`Area = Amp × Sigma × √(2π) / binWidth`), with the uncertainty propagated from the fit's full covariance matrix (including the Amp–Sigma correlation).

## Customizing a fit

All fit inputs are set as plain variables near the top of the macro's main function:

```cpp
const int    N_PEAKS     = 5;                            // number of peaks in the ROI
const bool   USE_EXP_BKG = false;                         // false = linear bkg, true = exponential
std::vector<double> PEAK_GUESS  = {879, 911, 940, 968, 1009};  // initial peak positions
std::vector<double> SIGMA_GUESS = {12, 5, 6, 5, 3.5};          // initial peak widths
const double BKG_LOW_MIN  = 820.0, BKG_LOW_MAX  = 855.0;  // peak-free window below the peaks
const double BKG_HIGH_MIN = 1035.0, BKG_HIGH_MAX = 1060.0; // peak-free window above the peaks
```

### Fixing individual peak parameters

Any parameter of any peak can be held fixed at a known value instead of being freely fitted, via the `FIXED_PARAMS` vector:

```cpp
FIXED_PARAMS = {
    { 3, ParamType::kMean,  940.0 },   // fix peak #3's position
    { 4, ParamType::kSigma, 5.0   },   // fix peak #4's width
};
```

## Notes

- Peak numbering in all output is 1-based (`Peak 1` = the first entry in `PEAK_GUESS`).
