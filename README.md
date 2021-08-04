# LyheeTools
A few rough solutions for astronomical data handling problems I met

## Directories
**TransientSurvey4yr**: Codes used for JCMT Transient Survey summary paper. The data used is confidential. Accepted for publication in *The Astrophysical Journal*. (https://arxiv.org/abs/2107.10750)

**EC53Paper**: Codes used for the analysis of EC53 using confidential data sets. Published in *The Astrophysical Journal*. Lee Y.-H. et al., 2020, ApJ, 903, 5 (https://ui.adsabs.harvard.edu/abs/2020ApJ...903....5L/abstract)


## Selfmade Packages

### LCAnalyses


- `AutoCorrelation(dates, fluxes, flux errors, (interval=30, custom=0))`
	- Resamples the lightcurve as given interval (or given window function) by linear interpolation and calculate the autocorrelation function (ACF) in given lag (set by given interval)
	- **Return**: resampling function, resampled fluxes, resampled noises, ACF
	- **Useful option**: You can put your own resampling function in the *custom* option

- `LCResample(dates, fluxes, flux errors, window function)`
	- Resamples the lightcurve as given window function by linear interpolation
	- **Basic Input**: dates, fluxes, flux errors, window function
	- **Return**: resampling dates, resampled fluxes, resampled noises

- `LCSEDflux(dates, fluxes, flux errors, resampling dates, wavelength)`
	- Calculates the SED flux in wavelength domain
	- **Return**: resampling dates, SED fluxes, SED flux errors

- `LCStringlength(dates, fluxes, flux errors, periods, (phaseplot))`
	- Measures the string-lengths of the phase diagrams obtained from given periods
	- **Return**: string lengths, string lengths error
	- **Useful option**: Giving a period in *phaseplot* will display the phase diagram with the string at a given period.

- `LCDiffCoeff(dates, fluxes, flux errors)`
	- Calculates the differential coefficients of the lightcurve between every adjacent data points
	- **Return**: dates in the middle, differential coefficients, errors

- `LCLombscargle(dates, fluxes, flux errors)`
	- Object for easier access to the Lomb-Scargle periodogram analysis
	- Can get Periodogram (frequency, statistical powers), best-fit sinusoid (or their parameters), and window periodogram

- `LCSampling(dates, fluxes, flux errors, rate=0.8)`
	- Randomly samples the given rate from the lightcurve (datapoints in general)
	- **Return**: sampled dates, sampled fluxes, sampled flux errors

- `LCLinearfit(data_x, data_y, error_y, (error_x))`
	- Finds the best-fit parameters of linear function for the dataset (using *scipy.optimize.curve_fit* in default)
	- **Return**: solution parameters, solution errors
	- **Useful option**: if *error_x* is given, it uses *scipy.odr*.

- `Gaussfit(data_x, data_y, error_y)`
	- Finds the best-fit parameters of Gaussian function for the dataset (using *scipy.optimize.curve_fit* in default)
	- **Return**: Gaussian parameters, parameter errors

- `LCPeriodogramerr(frequencies, powers, peakposition=0, plot=0)`
	- Performs the Guassian fitting on the peak of periodogram
	- **Return**: Gaussian parameters, parameter errors, error range of the peak period
	- **Useful option**: Default run finds the best-fit peak (highest peak)

- `LCTimeWindCut(dates, fluxes, flux errors, initial date, final date)`
	- Cuts the lightcurve following the given dates
	- **Return**: cut dates, cut fluxes, cut flux errors


*JDUTPlot.py*: Overplots UT year on top of the brightness-JD lightcurve, with the year separation grid.

*FITS2ASCII.py*: Make FITS format image data to ASCII data with three columns (RA[Deg], Dec[Deg], and Flux).
