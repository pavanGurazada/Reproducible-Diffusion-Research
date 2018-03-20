using LightGraphs, Distributions
using DataFrames
using ProgrssMeter

srand(20130810)

"""

In this script, we explore the results from "Temporal profiles of avalanches
on networks", Gleeson and Durrett (2017). We replicate the results based on the
Centola-Macy model as discussed in the paper. The central idea is that the
temporal distribution of a large number of cascades run for a fixed time T
can be expressed as ``f(t/T)``. While the paper derives closed-form solutions
for specific cases, in this script we wish to explore two specific questions:
    1. How long did each of these cascades last, i.e., what is the statistical
       distribution of ``f(t/T)``?
    2. How many nodes get engaged in these cascades?

We consider N = 10^3 avalanches, with the duration of avalanches considered to be
a maximum of T = 20.

Each node has an integer threshold, drawn at random from a uniform distribution

Every cascade is assumed to be initiated by a randomly-chosen single node, which
is activated at the beginning of the process.

Further to activation, cascades proceed according to the threshold model as
follows:

1. Each node i has a positive threshold R_i that is assigned randomly from a
distribution.

2. When an inactive node i of degree k_i is chosen for updating, it considers
the number m_i of its neighbors who are active, and if m_i >= R_i, it activates

For threshold models thetamax (maximum threshold for a node) = <k^2 - k>/z
"""

struct CentolaMacyModel
    G::LightGraphs.SimpleGraphs.SimpleGraph{Int}
    threshold::Vector{Int}
    node_status::BitVector
end

# For those networks close to criticality
function initialize_threshold!(M::CentolaMacyModel, theta_max::Float64) # theta_max = 5.0736
    append!(M.threshold, ceil.(Int, thetamax .* rand(nv(M.G))))
end

# logical definition for normal networks where the threshold can be any number
# from 1 to the number of neighbors of each node
function initialize_threshold!(M::CentolaMacyModel)

    for node in vertices(M.G)
        push!(M.threshold, sample(1:length(neighbors(M.G, node))))
    end

end

function seed!(M::CentolaMacyModel)
    seed = sample(vertices(M.G))
    M.node_status[seed] = true
end

function n_active_nbrs(M::CentolaMacyModel, node::Int)
    return sum([M.node_status[nbr] for nbr in neighbors(M.G, node)])
end

# Update the vulnerable nodes in the network using the their thresholds
function evolve!(M::CentolaMacyModel)
    for node in shuffle(vertices(M.G))
        if (M.node_status[node] == false) &&
           (n_active_nbrs(M, node) >= M.threshold[node])
           M.node_status[node] = true
       end
   end
end

# Run the simulations on the Barabasi-Albert networks
# How long did each of these cascades last?
# How many nodes get engaged in these cascades?

function simulate(n::Int, k::Int, n_realizations::Int; T = 20)

    end_times = DataFrame(r = Int[], end_time = Int[], engaged_nodes = Int[])

    @showprogress 1 "Simulating..." for r in 1:n_realizations

        g = barabasi_albert(n, k)

        M = CentolaMacyModel(g, Int[], falses(nv(g)))
        initialize_threshold!(M)
        seed!(M)
        evolve!(M)
        n_newly_active = sum(M.node_status)

        t = 0

        while (n_newly_active > 0) && (t < T)
            t += 1
            evolve!(M)
            n_newly_active = sum(M.node_status) - n_newly_active
        end

        push!(end_times, [r, t, sum(M.node_status)])

    end

    return end_times

end

function simulate(n::Int, k::Int, theta_max::Float64, n_realizations::Int; T = 20)

    end_times = DataFrame(r = Int[], end_time = Int[], engaged_nodes = Int[])

    @showprogress "Simulating..." for r in 1:n_realizations

        g = barabasi_albert(n, k)

        M = CentolaMacyModel(g, Int[], falses(nv(g)))
        initialize_threshold!(M, theta_max)
        seed!(M)
        evolve!(M)
        n_newly_active = sum(M.node_status)

        t = 0

        while (n_newly_active > 0) && (t < T)
            t += 1
            evolve!(M)
            n_newly_active = sum(M.node_status) - n_newly_active
        end

        push!(end_times, [r, t, sum(M.node_status)])
    end

    return end_times
end

@time results = simulate(10^5, 3, 100)
@time results_thetamax = simulate(10^5, 3, 5.0736, 1000)
