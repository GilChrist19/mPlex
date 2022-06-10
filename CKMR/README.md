# CKMR: Agent-based Modeling of Close-Kin Mark Recapture Methods in Mosquitoes
Building on the basic mPlex organization, this branch allows for sampling of individuals at every life stage to perform [close-kin analysis](https://projecteuclid.org/download/pdfview_1/euclid.ss/1464105042).  

## Software
  * inst/include/ contains header files
  * man contains auto-generated function documentation
  * R/ contains R files
  * src/ contains implementation files and makefiles
  * TESTFILE.R is a basic example of how to run this package

## Installation

There are several methods for installing packages. Because this package includes 
compiled code, there are additional dependencies that must be installed outside of 
R. 

### External Dependencies
  * Windows: The [Rtools](https://cran.r-project.org/bin/windows/Rtools/) toolchain 
    must be installed in order to build/compile packages. ([Additional Info](https://cran.r-project.org/bin/windows/base/howto-R-devel.html))
  * Mac: This packages requires `RcppArmadillo`, which has Fortran dependencies. 
    Macs have buggered their compiler toolchain, and it needs rebuilt to support 
    Fortran linking. Instructions to do this are [here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).
  * Linux: You should know what's going on, or you'll learn very quickly. 

### Package Dependencies

This package depends on several other packages to run ([Rcpp](https://www.rcpp.org/), 
[RcppProgress](https://cran.r-project.org/package=RcppProgress), [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)) and to build 
documentation ([knitr](https://yihui.org/knitr/), [rmarkdown](https://rmarkdown.rstudio.com/), 
[markdown](https://cran.r-project.org/package=markdown)).  
Additionally, the [devtools](https://devtools.r-lib.org/) library makes several 
installation steps significantly easier.  

These can be installed from within R using the following command:
`install.packages(pkgs = c('devtools','Rcpp','RcppArmadillo','RcppProgress','knitr','rmarkdown','markdown'))`

### Installation Methods
* Github Method
  * This is the easiest method, just requires `devtools` and `R`.
  1. `devtools::install_github(repo = "GilChrist19/mPlex/CKMR", build_manual = TRUE, 
  build_vignettes = TRUE, upgrade = "ask")`
* RStudio Method
  * This is the next easiest method, and can be done all in the GUI.
  1. Download this repository.
  2. Navigate into the package directory - this is the `~/MGDrivE2/` directory.
  3. Double click the `MGDrivE2.Rproj` file - this opens the packge in RStudio.
  4. Go the to `Build` tab at the top of the window.
  5. click 'load all'
  6. Click `clean and rebuild` - this will build and install the package, and the 
    basic documentation, but not the vignettes.
* CMDline Method (in R)
  * This is very similar to the RStudio method, and makes use of the `devtools` package.
  1. Download this repository.
  2. In the CMDline, navigate to the package directory  - this is the `~/MGDrivE2/` directory.
  3. Begin R, and type `devtools::document()` - this builds the basic documentation.
  4. When that is finished, type `devtools::build()` - this actually builds the package.
  5. Finaly, type `devtools::install()` - this finishes the installation process, 
    everything is ready to go.
* CMDline Method (bash)
  * This is a more difficult method, but does not involve extra packages such as `devtools`
  1. Download this repository.
  2. Navigate into the package directory - this is the `~/MGDrivE2/` directory.
  3. Type `R CMD build ./` - this builds the package and all documentation, including 
    vignettes by default.
  4. Type `cd ../` (alternatively, you can stay 1 level higher in step 2, and replace 
    `./` with `MGDrivE2`)
  5. Type `R CMD INSTALL MGDrivE2.tar.gz` - this finishes the installation.

### References/Additional Info
  * [MIT](https://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf)
  * [Coding Club](https://ourcodingclub.github.io/tutorials/writing-r-package/)
  * [Prestevez](https://www.prestevez.com/post/r-package-tutorial/)
  * [kbroman](https://kbroman.org/pkg_primer/pages/build.html)
  * I like everything this guy has blogged.
  * [r-pkgs](https://r-pkgs.org/index.html)
  * This one I can vouch for, it's also an entire book. 
  * [adv-r](https://adv-r.hadley.nz/)
  * This will take you zero-to-hero, but it's a deep dive.
  * [r4ds](https://r4ds.had.co.nz/index.html)
  * R for data science, also includes projects and packages.
  * [More Advanced](https://support.rstudio.com/hc/en-us/articles/200486488-Developing-Packages-with-the-RStudio-IDE)
  * Several links to other resources, not the most clear, but very complete.
  * Includes some info for windows vs linux.
  * [Video](https://www.youtube.com/watch?v=79s3z0gIuFU)
  * 15-min walkthrough of setting up a new package through RStudio.
  * [Slide Deck](https://johnmuschelli.com/smi_2019/index.html#1)
  
## Authors
Jared Bennett <jared_bennett@berkeley.edu>

## To-do


