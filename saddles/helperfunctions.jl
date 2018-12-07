using LightGraphs
using Random: shuffle, seed!

seed!(20130810)


function prob_activation(early_market_graph::LightGraphs.SimpleGraphs.SimpleGraph, 
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
        println(n_active_nbrs)
        n_active_early_nbrs = sum([early_node_status[early_nbr] 
                                   for early_nbr in rand(vertices(early_market_graph), n_weak_ties)])
        println(n_active_early_nbrs)
        return 1 - (1-probs["p_m"]) * ((1 - probs["q_mm"])^n_active_nbrs) * 
                                      ((1 - probs["q_im"])^n_active_early_nbrs)
    end
end

early_market_graph = erdos_renyi(1000, 0.5)
main_market_graph = erdos_renyi(1000, 0.2)

early_node_status = falses(nv(early_market_graph)) 
main_node_status = falses(nv(main_market_graph))

probs = Dict("p_i" => .01,
             "p_m" => .0005,
             "q_ii" => .005,
             "q_mm" => .0005,
             "q_im" => 5e-7)

println(prob_activation(early_market_graph, main_market_graph, 1, 5, 
                        early_node_status, main_node_status, probs, false))
