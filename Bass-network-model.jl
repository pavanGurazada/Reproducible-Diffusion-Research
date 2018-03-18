using LightGraphs, Distributions # To run the simulation
using DataFrames, RCall, JuliaDB # To analyze results

srand(20130810)

"""

In this script, we explore the modeling guidelines presented in Rand and Rust (2011). 

For the original Bass model, we use a complete graph assuming that each node is connected to all the other nodes. Hence, the behavior of all other users in a market is visible to each individual user who make their decisions based on this knowledge.

We then incorporate the impact of network structures by restricting the visibility of the nodes to their neighborhoods. We consider two types of graphs for this: Erdos-Renyi random networks and Barabasi-Albert scale free networks.

All simulations follow the recipe:

initialize(...) |>  evolve!(...) 

"""

struct BassModel
    G::LightGraphs.SimpleGraphs.SimpleGraph{Int}
    p::Float64
    q::Float64
end

function adoption_prob(node_status::BitVector,
                       M::BassModel, 
                       node::Int)

    n_adopted_nbrs = sum([node_status[nbr] for nbr in neighbors(M.G, node)])

    return M.q * n_adopted_nbrs/nv(M.G)

end

function initalize(M::BassModel)
    return falses(nv(M.G))
end

function evolve!(node_status::BitVector, M::BassModel)

    for node in shuffle(vertices(M.G))
        if node_status[node] == false
            if (rand(Uniform()) < M.p) || 
                (rand(Uniform()) < adoption_prob(node_status, M, node))

                node_status[node] = true
            end
        end
    end

    return nothing
end

# Run the simulations on the original Bass model

function simulate(M::Vector{Int}, P::Vector{Float64}, Q::Vector{Float64},                         n_realizations::Int)
    
    parameter_space = DataFrame(m = M, p = P, q = Q)

    output = DataFrame(m = Int[], p = Float64[], q = Float64[],
                       t = Int[], n_adopted = Int[])

    for params in eachrow(parameter_space)
        for r in n_realizations 
            
        end
    
    end


end