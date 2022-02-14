# Generalized Additive Models; a data-driven approach to estimating regression models

### Physalia-Courses 

https://www.physalia-courses.org/

### Gavin Simpson

#### 14&ndash;18th February, 2022

## Overview

Most of the statistical methods you are likely to have encountered will have specified fixed functional forms for the relationships between covariates and the response, either implicitly or explicitly. These might be linear effects or involve polynomials, such as x + x<sup>2</sup> + x<sup>3</sup>. Generalised additive models (GAMs) are different; they build upon the generalised linear model by allowing the shapes of the relationships between response and covariates to be learned from the data using splines. Modern GAMs, it turns out, are a very general framework for data analysis, encompassing many models as special cases, including GLMs and GLMMs, and the variety of types of splines available to users allows GAMs to be used in a surprisingly large number of situations. In this course we’ll show you how to leverage the power and flexibility of splines to go beyond parametric modelling techniques like GLMs.

## Target audience and assumed background

The course is aimed at at graduate students and researchers with limited statistical knowledge; ideally you’d know something about generalised linear models. But we’ll recap what GLMs are so if you’re a little rusty or not everything mentioned in the GLM course makes sense, we have you covered.

Participants should be familiar with RStudio and have some fluency in programming R code, including being able to import, manipulate (e.g. modify variables) and visualise data. There will be a mix of lectures, in-class discussion, and hands-on practical exercises along the course.

## Learning outcomes

 1. Understand how GAMs work from a practical view point to learn relationships between covariates and response from the data
 2. Be able to fit GAMs in R using the mgcv and brms packages
 3. Know the differences between the types of splines and when to use them in your models
 4. Know how to visualise fitted GAMs and to check the assumptions of the model

## Programme

Sessions from 14:00 to 20:00 (Monday to Thursday), 14:00 to 19:00 on Friday (Berlin time). From Tuesday to Friday, the first hour will be dedicated to Q&A and working through practical exercises or students’ own analyses over Slack and Zoom. Sessions will interweave mix lectures, in-class discussion/ Q&A, and practical exercises.

### Monday

[Slides](https://gavinsimpson.github.io/physalia-gam-course/index.html)

* Brief overview of R and the Tidyverse packages we’ll encounter throughout the course
* Recap generalised linear models
* Fitting your first GAM

### Tuesday

* How do GAMs work?
* What are splines?
* How do GAMs learn from data without overfitting?

We’ll dig under the hood a bit to understand how GAMs work at a practical level and how to use the mgcv and gratia packages to estimate GAMs and visualise them.

### Wednesday

* Model checking, selection, and visualisation.
* How do we do inference with GAMs?
* Go beyond simple GAMs to include smooth interactions and models with multiples smooths.

### Thursday

* Hierarchical GAMs; introducing random smooths and how to model data with both group and individual smooth effects.
* Doing more with your models; introducing posterior simulation.

### Friday

* Fitting GAMs with brms
* Going beyond the mean; fitting distributional models and quantile GAMs

