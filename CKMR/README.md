# CKMR: Agent-based Modeling of Close-Kin Mark Recapture Methods in Mosquitoes
Building on the basic mPlex organization, this branch allows for sampling of individuals at every life stage to perform [close-kin analysis](https://projecteuclid.org/download/pdfview_1/euclid.ss/1464105042).  

## Software
  * inst/include/ contains header files
  * src/ contains implementation files and makefiles
  * R/ contains R files
  * man/ contains automatically generated R documentation files

## Authors
Jared Bennett, [Sean Wu](https://slwu89.github.io), [Héctor Manuel Sánchez Castellanos](https://chipdelmal.github.io), and [John M. Marshall](http://sph.berkeley.edu/john-marshall)

## To-do
  1. Write out pseudocode of `offspringDistribution` function.
    * goal is to de-couple the tightly coupled `reference` and `offspringDistribution` elements; we want a generic `offspringDistribution` so the only variant type is `reference`.
  2. Check move semantics between stl containers working properly for `Mosquito` for low overhead migration functions.
