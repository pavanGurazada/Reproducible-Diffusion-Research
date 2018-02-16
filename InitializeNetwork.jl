using LightGraphs, SimpleWeightedGraphs
using Distributions

srand(20130810)

function count_ties(G::SimpleWeightedGraphs.SimpleWeightedGraph, node::Int)
    str_ties, wk_ties = 0, 0

    for nbr in neighbors(G, node)
        wt = get_weight(G, node, nbr)

        if wt == 0x2
            str_ties += 1
        end

        if wt == 0x1
            wk_ties += 1
        end
    end
    return (str_ties, wk_ties)
end

function initialize_network(n_nodes::Int, n_strong_ties::Int, n_weak_ties::Int)
    G = SimpleWeightedGraph(random_regular_graph(n_nodes, n_strong_ties + n_weak_ties))

    for node in vertices(G)
        str_ties, wk_ties = count_ties(G, node)

        if str_ties < n_strong_ties
            str_nbrs = sample(neighbors(G, node),
                              n_strong_ties - str_ties,
                              replace = false)
            for nbr in str_nbrs
                add_edge!(G, node, nbr, 0x2)
            end
        end

        str_ties, wk_ties = count_ties(G, node)

        if wk_ties < n_weak_ties
            wk_nbrs = sample(neighbors(G, node),
                             n_weak_ties - wk_ties,
                             replace = false)
            for nbr in wk_nbrs
                add_edge!(G, node, nbr, 0x1)
            end
        end
    end

    return G
end

@time g = initialize_network(3000, 25, 5)

str_degrees, wk_degrees = Int[], Int[]
for node in vertices(g)
    push!(str_degrees, count_ties(g, node)[1])
    push!(wk_degrees, count_ties(g, node)[2])
end

mean(str_degrees)
mean(wk_degrees)
