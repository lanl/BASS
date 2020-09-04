# BASS 1.2.1
- remove dependency on depricated wtmsa package and remove wavelet example
- add reference to JSS paper

# BASS 1.2.1
- change dependency from gsl to hypergeo so that the user does not need gsl

# BASS 1.2.0
- Sobol getEffects option added
- changed how Wallenius' noncentral hypergeometric density is calculated (in R rather than C++ from BiasedUrn)
- various bug fixes

# BASS 1.1.2
- added unit tests for prediction and Sobol decomposition for various validated models
- added ability to write models to JSON, and read models from JSON
- various options added to existing functions
- existing function bassOB renamed to bassBasis

# BASS 1.1.1
- added non-uniform prior support for continuous variables in Sobol decomposition

# BASS 1.1.0
- fixed bug in acceptance ratios
- added support for BASS in principal component space, including Sobol decomposition

# BASS 0.2.3 
- added a lower bound parameter for the prior on s2 (the error variance)

# BASS 0.2.2 
- added print and summary methods for bass objects

# BASS 0.2.1
- added "small" parameter to allow for smaller memory footprint
- removed uses of timestamp to avoid a bug in Rstudio on Windows
- corrected example to use temperatures rather than inverse temperatures, and to match vignette

# BASS 0.2.0
- vignette added
- argument changes
    - temp.ladder now refers to temperatures rather than inverse temperatures
- Bug fixes:
    - Sobol index rounding for functional case
    - Handling responses in a dataframe
    - Error handling for w1 or w1 equal to 0
- Documentation updates
    - Changed input names a.beta.prec and b.beta.prec from the g-prior to a.tau and b.tau

# BASS 0.1.0
- First release.
