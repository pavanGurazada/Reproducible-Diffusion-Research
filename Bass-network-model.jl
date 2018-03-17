using LightGraphs, Distributions # To run the simulation
using DataFrames, RCall, JuliaDB

srand(20130810)

"""

In this script, we explore the modeling presented in Rand and Rust (2011). 
The excellent guidelines presented in this paper are implemented in this paper
"""

struct Network
    G::LightGraphs.SimpleGraphs.SimpleGraph{Int}
    p::Float64
    q::Float64
end

