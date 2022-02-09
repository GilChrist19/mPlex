# mPlex: Agent-based Modeling of Multiplex Gene Drive in Mosquitoes
mPlexCpp is a C++11/14-based R package to run simulations of multiplex (including daisy-drive) gene drive systems. It is being developed at the [Marshall Lab](https://www.marshalllab.com) at the University of California, Berkeley.

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

## Notes and References
 1. Short and sweet [OMP guide](https://chryswoods.com/beginning_openmp/README.html)
