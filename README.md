# Generalized Additive Models; a data-driven approach to estimating regression models

### Physalia-Courses 

https://www.physalia-courses.org/

### Gavin Simpson

#### 8&ndash;11th December, 2025

## Overview

Most of the statistical methods you are likely to have encountered will have specified fixed functional forms for the relationships between covariates and the response, either implicitly or explicitly. These might be linear effects or involve polynomials, such as x + x<sup>2</sup> + x<sup>3</sup>. Generalized additive models (GAMs) are different; they build upon the generalized linear model (GLM) by allowing the shapes of the relationships between response and covariates to be learned from the data using splines. Modern GAMs, it turns out, are a very general framework for data analysis, encompassing many models as special cases, including GLMs and GLMMs, and the variety of types of splines available to users allows GAMs to be used in a surprisingly large number of situations. In this course we’ll show you how to leverage the power and flexibility of splines to go beyond parametric modelling techniques like GLMs.

## Target audience and assumed background
The course is aimed at at graduate students and researchers with limited statistical knowledge; ideally you’d know something about generalized linear models, but we’ll recap what GLMs are, so if you’re a little rusty or not everything mentioned in a GLM course made sense, we have you covered.

Participants should be familiar with RStudio and have some fluency in programming R code, including being able to import, manipulate (e.g. modify variables) and visualise data. There will be a mix of lectures, in-class discussion, and hands-on practical exercises along the course. From running the course previously, knowing the difference between "fixed" and "random" effects, and what the terms "random intercepts" and "random slopes" are, will be helpful for the Hierarchical GAM topic, but we don't expect you to be an expert in mixed effects or hierarchical models to take this course.

## Learning outcomes

1. Understand how GAMs work from a practical view point to learn relationships between covariates and response from the data,

2. Be able to fit GAMs in R using the mgcv package,

3. Know the differences between the types of splines and when to use them in your models,

4. Know how to visualise fitted GAMs and to check the assumptions of the model.

## Pre-course preparation

### Install an up-to-date version of R

Please be sure to have at least version 4.5.0 of R installed (the version of my gratia package we will be using depends on you having at least version 4.1.0 installed and some slides might contain code that requires version 4.4.x). Note that R and RStudio are two different things: it is not sufficient to just update RStudio, you also need to update R by installing new versions as they are release.

To download R go to the [CRAN Download](https://cran.r-project.org/) page and follow the links to download R for your operating system:

* [Windows](https://cran.r-project.org/bin/windows/)
* [MacOS X](https://cran.r-project.org/bin/macosx/)
* [Linux](https://cran.r-project.org/bin/linux/)

To check what version of R you have installed, from within R, you can run

```r
version
```

then look at the `version.string` entry (or the `major` and `minor` entries). For example, on my system I see:

```
# ... output not shown ...
major          4                           
minor          5.1
# ... output not shown ...
version.string R version 4.5.1 (2025-06-13)
# ... output not shown ...
```

### Update your R packages, and install the required R packages

We will make use of several R packages that you'll need to have installed. Prior to the start of the course, please run the following code to update your installed packages and then install the required packages:

```r
# update any installed R packages
update.packages(ask = FALSE, checkBuilt = TRUE)

# packages to install
pkgs <- c("gamm4", "tidyverse", "readxl", "mgcViz", "DHARMa", "gratia",
  "marginaleffects", "ggforce")

# install those packages
install.packages(pkgs, Ncpus = 4) # set Ncpus to # of *physical* CPU cores you have
```

We might need the development version of *gratia*; you can install this with

```r
# Install gratia in R
install.packages("gratia", repos = c(
  "https://gavinsimpson.r-universe.dev",
  "https://cloud.r-project.org"
))
```

Now we must check that we actually do have recent versions of the packages installed; if your R is not reasonably new (gratia requires R>= 4.1.0, but some of the *tidyverse* packages may need an R that is newer than this) then you may be stuck on out-dated versions of the packages listed above. This is why I recommend that you install the latest version of R. If you choose to use an older version of R than version 4.5.x (where *x* is 0, 1, or 2 currently) then you do so at your own risk and you cannot expect support with setup problems during the course.

```r
vapply(c("mgcv", pkgs), packageDescription, character(1), drop = TRUE,
  fields = "Version")
```

On my system I see:

```r
> vapply(c("mgcv", pkgs), packageDescription, character(1), drop = TRUE,
    fields = "Version")
         mgcv           gamm4       tidyverse          readxl          mgcViz
      "1.9-4"         "0.2-7"         "2.0.0"         "1.4.5"         "0.2.1"
       DHARMa          gratia marginaleffects
      "0.4.7"        "0.11.1"        "0.31.0"
```

The key ones are to be sure that *gratia* is version "0.11.1", *mgcv* is at least "1.9-4", and *tidyverse* is "2.0.0".

<!-- ### Installing the *cmndstan* backend (optional)

Fitting GAMs with Stan is quite time consuming if we use the standard *rstan* interface. To speed things up significantly, we can use the *cmdstan* backend, however this requires a little more setup. If you can't get this to work don't worry, it's not an integral part of the course, as you can still use the *rstan* backend with `brm()`.

*cmdstan* requires a working C++ compiler on your system. Typically, Windows and MacOS X machines do not come with one installed by default. To install the C++ toolchain required you should follow the instructions [here](https://mc-stan.org/docs/cmdstan-guide/installation.html#cpp-toolchain), only the bits in the **C++ Toolchain** section that is linked to. If you're on a recent MacOS X system, installation of the required toolchain is relatively simple, requiring only installation of some parts of *xcode*. On Windows, things are slightly more complicated as you need to install the version of RTools for your version of R and then add some details to your `PATH` to allow the toolchain to be run from the command line. There are slightly different instructions (versions of RTools) to install depending on your version of R. *If any of this sounds too complicated for you, just stop here and don't proceed; you don't need to run the `brm()` code when I am working through some examples and we won't spend a lot of time on fully Bayesian GAMs anyway.*

Once you have the toolchain installed, to do the actual installation of the *cmdstan* backend we need to load the *cmndstanr* package and complete some steps. Give yourself some time to do this as the options below will download the backend and start to compile it for your computer.

```r
# install cmdstanr
install.packages("cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# load the R package interface to cmdstan
library(cmdstanr)

# check the your toolchain is configured correctly and working
check_cmdstan_toolchain()
# if this says anything other than that the toolchain is configured properly
# stop(!) and go back to the C++ Toolchain instructions and make sure you
# have completed all the steps for your OS

# install cmndstan backend
# You can increase `cores` if you have more cores available on your system
# if in doubt, just leave it as shown below
install_cmdstan(cores = 2)
# wait for some time...

# you can confirm that cmndstan is installed and what version you have with
cmdstan_path()
cmdstan_version()
```
-->

<!-- Finally, we will make use of the development version of the gratia package as it is not quite ready for CRAN. You can install this package using the binaries provided by the [rOpenSci](https://ropensci.org/) build service [R-Universe](https://r-universe.dev). To install from my R-Universe, you need  to tell R to also install packages from my R-Universe package repo:

```r
# Download and install gratia
install.packages("gratia",
    repos = c("https://gavinsimpson.r-universe.dev", "https://cloud.r-project.org"))

```
-->

## Programme

Sessions from 14:00 to 20:00 (Monday to Thursday), Berlin time. Sessions will interweave mix lectures, in-class discussion/ Q&A, and practical exercises.

### Monday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/day-1/index.html)

* Brief overview of R and the Tidyverse packages we’ll encounter throughout the course
* Recap generalised linear models
* Fitting your first GAM

### Tuesday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/day-2/index.html)

* How do GAMs work?
* What are splines?
* How do GAMs learn from data without overfitting?

We’ll dig under the hood a bit to understand how GAMs work at a practical level and how to use the mgcv and gratia packages to estimate GAMs and visualise them.

### Wednesday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/day-3/index.html)

* Model checking, selection, and visualisation.
* How do we do inference with GAMs?
* Go beyond simple GAMs to include smooth interactions and models with multiples smooths.

### Thursday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/day-4/index.html)

* Hierarchical GAMs; introducing random smooths and how to model data with both group and individual smooth effects.
* Doing more with your models; introducing posterior simulation.
* Worked examples

<!-- ### Friday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/day-5/index.html)

* Going beyond the mean; fitting distributional models

* Worked examples -->
