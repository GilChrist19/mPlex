# Verily Planning Simulation

Verily has collected mosquitoes for several prior studies, but needs input on 
what mosquitoes to sequence for our project. The idea is, we suggest a small sample 
of mosquitoes to start with, they sequence them, we analyze them, then based on that 
sample and analysis we request many more mosquitoes sequenced, Verily sequences 
the second set, and then we do a final analysis using all of the samples.  
  
This folder sets up the initial simulations to determine what mosquitoes to start 
with. It will include several files:

* [Locations](./c2_centroids_info.rds)
  * This file contains lat/longs for all of the houses and traps in Verily's C2
  control area, from the [Fresno Study](https://doi.org/10.1038/s41587-020-0471-x).  
  The file has 4 columns: ID, longitude, latitude, trap (boolean).
* [Kernel1](./c2_kernel1.rds)
  * This is a kernel that Tomas generated, that corresponds to the trapping numbers 
  from the Fresno Study.
* [Kernel2](./c2_kernel_exp80.rds)
  * Synthetic kernel, assuming an 80% life-time probability to stay in your original 
  location, an average movement distance (given that you do move) of 112.5m, and 
  an adult mortality of 0.09 per day.
* [kernelGen.R](./kernelGen.R)
  * This file generates kernels from the [locations](./c2_centroids_info.rds) file. 
  This way, we can create specific kernels, depending on what we want to test. 
  These will be less realistic than the kernel from Tomas, but more uniform.
* [Auxiliary File 1](./combineFiles.R)
  * This file supplies 2 functions for combining the output from CKMR simulations. 
  It should not be modified here, as the original design is elsewhere. It is only 
  included for compeleteness.
* [Simulations](./sims.R)
  * This file runs, cleans, and stores simulations over different kernels. The 
  output is stored by date, then by kernel name, as tar.bz2 files. It manages all 
  directories required. Requires all of the kernel files and the auxiliary file.

### 20200902

Initial setup of this file, and the directory, for generating  simulated data to 
inform the first round of Verily sequencing. The files at the top of this script 
will be updated as things are added, notes below the dates will not be modified 
after the day.  
  
Moved locations and kernel files from Tomas here. Compressed them for faster loading 
in R and smaller size.  
  
Started the kernel generation script. It generates a plot of the locations, highlighting 
traps in the map. It also generates the exponential kernel for testing, based on 
distances, desired migration, and current life parameters.

### 20200903

Finished generating the first exponential kernel. Moving auxiliary files here 
for completeness.  
  
Created and tested the simulation script. It runs, properly manages files, and 
returns compressed data. The data is labeled by date, then kernel name. Need to 
check life and sampling params, then we're good to go.  
  
Update simTime to 190 days, 100 days of burn-in and 90 days, approximately 3 months, 
of data.  
All other parameters are good. Need to run over 2 adult mortalities (done).  
Data from today is in the [google drive](https://drive.google.com/drive/folders/1VVGG048C4giHtuDdN8bSyKAlUSm0-7Gw?usp=sharing).








