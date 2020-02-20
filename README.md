[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

# surfspec: UGent PAE tools for imaging spectroscopy of macrophytes and sediments

This package contains a bundle of R functions developed for internal use at the Protistology and Aquatic Ecology (PAE) research group from Ghent University, Belgium. Its development is not intended for public release of an R package, but is mainatined in GitHub for collaboration and documentation of methods. The functions are developed for improved calculation of in air Lambertian equivalent bi-hemispherical reflectance of macrophytes and sediment cores, with specific functions for sediment reflectance profiles. Current instrumentation at PAE consist of an imaging spectrometer (SOC-710P, Surface Optics Corporation), a NIST traceable grey card and natural illumination.

### Install from Github:

This library depends on packages sp, raster, rgdal, magrittr and igraph.

```
# install.packages(c("remotes", "sp", "raster", "rgdal", "magrittr", "igraph"))
remotes::install_github("AlexCast/surfspec")
```


