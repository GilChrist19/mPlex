# Agent-based Modeling of Mosquitoes

This repository contains 2 models that were developed at the [Marshall Lab](https://www.marshalllab.com) 
at the University of California, Berkeley.  
  
  * mPlex - an agent-based model for multiplexing (one locus, multiple gRNA or 
  multiple locus, single gRNAs at each) gene drives and daisy-drives. 
    * [mPlexR](./mPlexR): Pure R-based implementation of mPlex  
    This was the original test bed and is now out-of-date
    * [mPlexCpp](./mPlexCpp): Current, C++11/14-based package
  * [CKMR](./CKMR) - an agent-based adaptation of mPlex to study close-kin mark-recapture 
  methods in mosquitoes. This adaptation removes much of the genetics and allows for 
  more appropriate population sampling and geneology tracing. 

## Repository
  * [CKMR](./CKMR) - C++ package for close-kin mark-recapture studies
  * [Main](./Main) - Project directories and scripts
  * [mPlexCpp](./mPlexCpp) - C++ package for multiplex studies
  * [mPlexR](./mPlexR) - R package for multiplex studies (old, consider this deprecated)
  * [ComparisonScript](./ComparisonScript.R) - Comparison of mPlex vs [MGDrivE](https://cran.r-project.org/package=MGDrivE) 
  to check population dynamics and simple one-locus homing constructs. 

## Authors
Jared Bennett, [Sean Wu](https://slwu89.github.io), [Héctor Manuel Sánchez Castellanos](https://chipdelmal.github.io), and [John M. Marshall](http://sph.berkeley.edu/john-marshall)

## To-do
  1. Write out pseudocode of `offspringDistribution` function.
      * goal is to de-couple the tightly coupled `reference` and `offspringDistribution` 
    elements; we want a generic `offspringDistribution` so the only variant type is `reference`.
  2. Check move semantics between stl containers working properly for `Mosquito` for low overhead migration functions.
  3. Should steal the output/input classes from **MGDrivECpp** for more performance.

## Notes and References
 1. Short and sweet [OMP guide](https://chryswoods.com/beginning_openmp/README.html)
 2. [RcppArmadillo](http://arma.sourceforge.net/docs.html#uword) reference
 3. [Rcpp](https://teuder.github.io/rcpp4everyone_en/) reference

