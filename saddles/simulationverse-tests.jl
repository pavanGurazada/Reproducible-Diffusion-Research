include("simulationverse.jl")

early_market_graph = erdos_renyi(1000, 0.5)
main_market_graph = erdos_renyi(1000, 0.2)

early_node_status = falses(nv(early_market_graph))
main_node_status = falses(nv(main_market_graph))

n_weak_ties = 5

prob_params = Dict("p_i" => .01,
                   "p_m" => .0005,
                   "q_ii" => .005,
                   "q_mm" => .0005,
                   "q_im" => 5e-7)


evolve!(early_market_graph, main_market_graph, n_weak_ties, early_node_status,
        main_node_status, prob_params)

sum(early_node_status)
sum(main_node_status)
