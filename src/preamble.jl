#= ==========================================================================================
=============================================================================================

preamble

The environment is activated, all packages are loaded, and target paths for 
figures and animations are created.

=============================================================================================
========================================================================================== =#

#---------------------------------------- environment --------------------------------------#

using Pkg
Pkg.activate("rcsEnvironment")

#----------------------------------------- graphics ----------------------------------------#

using CairoMakie

#-------------------------------------- file management ------------------------------------#

using Dates
(!isdir("animations")) ? (mkdir("animations")) : nothing; (!isdir("animations/$(today())")) ? (mkdir("animations/$(today())")) : nothing;
(!isdir("figures")) ? (mkdir("figures")) : nothing; (!isdir("figures/$(today())")) ? (mkdir("figures/$(today())")) : nothing;