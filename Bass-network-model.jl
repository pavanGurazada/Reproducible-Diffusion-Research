using LightGraphs, Distributions # To run the simulation
using DataFrames, Query, RCall # To analyze results

srand(20130810)

"""

In this script, we explore the modeling guidelines presented in Rand and Rust (2011). 

For the original Bass model, we use a complete graph assuming that each node is 
connected to all the other nodes. Hence, the behavior of all other users in a market 
is visible to each individual user who make their decisions based on this knowledge.

We then incorporate the impact of network structures by restricting the visibility 
of the nodes to their neighborhoods. 

We consider two types of graphs for this: Erdos-Renyi random networks and 
Barabasi-Albert scale free networks.

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

function initialize(M::BassModel)
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

function simulate(Industry::Vector{String}, M::Vector{Int}, P::Vector{Float64},
                  Q::Vector{Float64}, n_realizations::Int)
    
    parameter_space = DataFrame(industry = Industry, m = M, p = P, q = Q)

    output = DataFrame(Industry = String[], m = Int[], p = Float64[], 
                       q = Float64[], r = Int[], t = Int[], n_adopted = Int[])

    for params in eachrow(parameter_space)
        for r in 1:n_realizations 
            M = BassModel(CompleteGraph(params[2]), params[3], params[4])
            node_status = initialize(M)
            for t = 1:30
                evolve!(node_status, M)
                push!(output, [params[1], params[2], params[3], params[4], r, t, 
                      sum(node_status)])
            end
        end
    end

    return output
end

const Industry = ["Refrigerators", "Freezers", "TV", "Softener", "AC", "Dryer", 
                  "Lawnmowers", "Bed", "Coffee", "Iron", "Player"]

const M = [1000 for i in 1:length(Industry)]

const P = [0.002617, 0.0181190, 0.0278770, 0.0177030, 0.0103990, 0.01712060, 
           0.0091837, 0.0058760, 0.0171350, 0.0286320, 0.0247960]

const Q = [0.21566, 0.17110, 0.251105, 0.29695, 0.41861, 0.35688, 0.33790, 
           0.24387, 0.30145, 0.32791, 0.65410]

@time results = simulate(Industry, M, P, Q, 10)

summary_results = by(results, 
                     [:Industry, :m, :p, :q, :t], 
                     df -> mean(df[:n_adopted]))
head(summary_results)

results_summary = @from i in results begin
                  @group i.n_adopted by [i.Industry, i.m, i.p, i.q, i.t] into g 
                  @select {mean_adopted = mean(g..n_adopted)}
                  @collect DataFrame
end