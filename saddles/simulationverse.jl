using LightGraphs
using Distributions

using Random: shuffle, seed!

seed!(20130810)

"""
`prob_activation` is a key function that computes the probability of activation of
a node based on its location, i.e., in the early market or the main market. A key
component is the execution of random meetings between the main market nodes and
the early market nodes. To execute this, we select a random sample from the nodes
in the early market and then check their status.
"""

function compute_activation_prob(early_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                                 main_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                                 node::Int,
                                 n_weak_ties::Int,
                                 early_node_status::BitVector,
                                 main_node_status::BitVector,
                                 probs::Dict,
                                 early_market_flag::Bool)


    if early_market_flag
        n_active_nbrs = sum([early_node_status[nbr] for nbr in neighbors(early_market_graph, node)])
        return 1 - (1-probs["p_i"]) * (1 - probs["q_ii"])^n_active_nbrs

    else
        n_active_nbrs = sum([main_node_status[nbr] for nbr in neighbors(main_market_graph, node)])
        n_active_early_nbrs = sum([early_node_status[early_nbr]
                                   for early_nbr in rand(vertices(early_market_graph), n_weak_ties)])
        return 1 - (1-probs["p_m"]) * ((1 - probs["q_mm"])^n_active_nbrs) *
                                      ((1 - probs["q_im"])^n_active_early_nbrs)
    end
end


"""
Evolve the status of all the nodes in the two networks by one step
"""

function evolve!(early_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                 main_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                 n_weak_ties::Int,
                 early_node_status::BitVector,
                 main_node_status::BitVector,
                 probs::Dict)

    for node in shuffle(vertices(early_market_graph))
        prob_active = compute_activation_prob(early_market_graph, main_market_graph,
                                              node, n_weak_ties, early_node_status, main_node_status,
                                              probs, true)

        if rand(Uniform()) < prob_active
            early_node_status[node] = true
        end
    end

    for node in shuffle(vertices(main_market_graph))
        prob_active = compute_activation_prob(early_market_graph, main_market_graph,
                                              node, n_weak_ties, early_node_status, main_node_status,
                                              probs, false)

        if rand(Uniform()) < prob_active
            main_node_status[node] = true
        end
    end

    return nothing
end
