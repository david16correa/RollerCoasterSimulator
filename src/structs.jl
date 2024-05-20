#= ==========================================================================================
=============================================================================================

structs

The structures needed for the simulations are defined

=============================================================================================
========================================================================================== =#

mutable struct curveStruct
    start::Float64
    stop::Float64
    length::Int64
    r::Function
end

mutable struct cartStruct
    s::Float64
    pos::Tuple{Float64, Float64, Float64}
    totalEnergy::Float64
    curve::curveStruct
end