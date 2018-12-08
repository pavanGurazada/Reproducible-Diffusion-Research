include("simulationverse.jl")


function simulate(early_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                  main_market_graph::LightGraphs.SimpleGraphs.SimpleGraph,
                  n_weak_ties::Int,
                  prob_params::Dict)

    early_node_status = falses(nv(early_market_graph))
    main_node_status = falses(nv(main_market_graph))

    total_engaged_nodes = sum(early_node_status) + sum(main_node_status)
    total_nodes = nv(early_market_graph) + nv(main_market_graph)

    engaged_nodes_timeseries = []

    while total_engaged_nodes <= .95 * total_nodes

        evolve!(early_market_graph, main_market_graph, n_weak_ties, early_node_status,
                main_node_status, prob_params)

        total_engaged_nodes = sum(early_node_status) + sum(main_node_status)

        append!(engaged_nodes_timeseries, total_engaged_nodes)

    end

    return engaged_nodes_timeseries

end
