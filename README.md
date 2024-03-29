# Repository for "Sacrum accelerations can predict whole-body kinetic and stride kinematics across running speeds"

Open Access Publication 👉 https://peerj.com/articles/11199/

## Project Contents
- **analysis.R** - This is probably the file you're looking for. It trains models, produces RMSE and
 MAPE values, and generates the subplots used in Figure 1. 
 
- **GRF_IMU_data.csv** - Data used to train/test Quantile Regression Forest and Linear Regression models.
    
    | Column Header | Contents (units) |
    | ------------- | ---------------- |
    | `Sub`         | Subject Number |
    | `IsFemale`    | Boolean indicating Subject Sex (0 = Male, 1 = Female) |
    | `Mass`        | Body Mass (kg) |
    | `Speed`       | Treadmill speed (m/s) |
    | `GRFPeak`     | Mean Peak vGRF calculated from treadmill data (body weight (BW)) |
    | `GRFImpulse`  | Mean Vertical Impulse calculated from treadmill data (BW*s) |
    | `GRFtc`       | Mean Ground Contact Time calculated from treadmill data (s) |
    | `IMUPeak`     | Mean Peak vGRF calculated from sacrum acceleration data (BW) |
    | `IMUImpulse`  | Mean Vertical Impulse calculated from sacrum acceleration data (BW*s) |
    | `IMUtc`       | Mean Ground Contact Time calculated from sacrum acceleration data (s) |
    | `IMUsf`       | Step frequency calculated from sacrum acceleration data (Hz) |
  
## Reproducibility Instructions
The analysis was performed with `R version 3.6.3 (2020-02-29)`. Reproducibility results may vary with different versions,
in particular those < `3.6.0`. See [this stack overflow post](https://stackoverflow.com/questions/47199415/is-set-seed-consistent-over-different-versions-of-r-and-ubuntu)
and the documentation for `sample.kind` and `RNGversion` [here](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html).

To replicate the environment used to perform the analysis, run the first couple of commented out lines in `analysis.R`
to install the version of the packages used to produce the results and figure in the manuscript:
```
# Install versioned packages used for analysis (if needed): ----
# A few sub-dependencies that don't always cooperate with install_version():
# install.packages('gower')
# install.packages('systemfonts')
# install.packages('gdtools')
# Primary dependencies:
# install.packages(remotes)
# library('remotes')
# install_version('caret', version = '6.0-85', dependencies = TRUE, repos = 'http://cran.us.r-project.org')
# install_version('quantregForest', version = '1.3-7', dependencies = TRUE, repos = 'http://cran.us.r-project.org')
# install_version('tidyr', version = '1.0.2', dependencies = TRUE, repos = 'http://cran.us.r-project.org')
# install_version('ggthemes', version = '4.2.0', dependencies = TRUE, repos = 'http://cran.us.r-project.org')  # may need {gdtools} and/or {systemfonts}
# install_version('ggridges', version = '0.5.2', dependencies = TRUE, repos = 'https://cran.us.r-project.org')
```

There may be additional dependencies that should be loaded via `install.packages()`. If prompted to upgrade packages, 
skip and press enter (this may take a while):
```
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

1: ALL
2: CRAN packages only
3: None
.
.
.
Enter one or more numbers, or an empty line to skip updates:
```
Running the rest of the script from source will produce the reported MAPE/RMSE values and the subplots in Figure 1.

## Contributions
- [Ryan Alcantara, University of Colorado](https://twitter.com/Ryan_Alcantara_)
- Evan Day, University of Oregon
- Mike Hahn, University of Oregon
- Alena Grabowski, University of Colorado
