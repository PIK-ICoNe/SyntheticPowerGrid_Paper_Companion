Pkg.activate(joinpath(@__DIR__, "../"))
using Revise
using SyntheticPowerGridsPaper

using PowerDynamics
using OrdinaryDiffEq
using PowerGridNoise  
using Interpolations
using Statistics
using DelimitedFiles

##
