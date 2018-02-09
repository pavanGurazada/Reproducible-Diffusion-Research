using Distributions
using LightGraphs

using Plots
plotlyjs()

srand(20130810)

function fraction_engaged(node::Int64,
                          G::LightGraphs.SimpleGraphs.SimpleGraph,
                          node_status::Vector{Int64})
    num_engaged_neighbors = 0
    for nbr in neighbors(G, node)
        if node_status[nbr] == 1
            num_engaged_neighbors += 1
        end
    end
    return num_engaged_neighbors/length(neighbors(G, node))
end

function update_node_status(G::LightGraphs.SimpleGraphs.SimpleGraph,
                            threshold::Vector{Float64},
                            node_status::Vector{Int64})

    new_node_status = copy(node_status)

    for node in vertices(G)
        if node_status[node] == 0
            if fraction_engaged(node, G, node_status) > threshold[node]
                new_node_status[node] = 1
            end
        end
    end
    return new_node_status
end

function update_node_status!(G::LightGraphs.SimpleGraphs.SimpleGraph,
                             threshold::Vector{Float64},
                             node_status::Vector{Int64})
    for node in vertices(G)
        if node_status[node] == 0
            if fraction_engaged(node, G, node_status) > threshold[node]
                node_status[node] = 1
            end
        end
    end

    return nothing
end

function diffusion_simulation(n::Int64,
                              k::Int64,
                              threshold::Vector{Float64},
                              T::Int64,
                              n_realizations::Int64,
                              mode::String)
    """
    This is the function we call to execute the diffusion simulation.
    It takes a graph, a vector of threshold fractions for each node,
    a initial seed fraction to initialize the network, and the number
    of time steps for which the diffusion simulation will be run

    The idea is to run the diffusion simulation a very large number
    of times and see in how many of these simulations we observe a
    global cascade, i.e., number of nodes engaged after the simulation
    process is a certain size of the network

    Hyper Parameters of the model
    ----------
    1. Threshold (distribution or a specific value)
    2. Number of nodes in the network (n)
    3. Edges attached each time (k)
    4. Synchronous or asynchronous updates
    """

    output = Vector{Int64}(n_realizations)

    if mode == "async"
        for r in 1:n_realizations
            G = barabasi_albert(n, k)
            # Select a single random node from the network and seed it
            node_status = zeros(Int64, nv(G))
            node_status[sample(vertices(G))] = 1

            for _ in 1:T
                update_node_status!(G, threshold, node_status)
            end
            output[r] = sum(node_status)
        end
    else
        for r in 1:n_realizations
            G = barabasi_albert(n, k)
            # Select a single random node from the network and seed it
            node_status = zeros(Int64, nv(G))
            node_status[sample(vertices(G))] = 1
            # Easier version
            for _ in 1:T
                node_status = update_node_status(G, threshold, node_status)
            end

            output[r] = sum(node_status)
        end
    end

    return output

end

function diffusion_simulation(n::Int64, # Number of nodes
                              k::Int64, # number of barabasi albert edges to rewire
                              threshold::Vector{Float64},
                              n_realizations::Int64)

    output = Vector{Int64}(n_realizations)

    for r in 1:n_realizations
        G = watts_strogatz(n, k, 0.5)
        # Select a single random node from the network and seed it
        node_status = zeros(Int64, nv(G))
        node_status[sample(vertices(G))] = 1

        new_node_status = update_node_status(G, threshold, node_status)
        t = 1

        # Keep updating node status till there are more nodes to activate
        while maximum(new_node_status .!= node_status) > 0
            node_status = new_node_status
            new_node_status = update_node_status(G, threshold, node_status)
            t += 1
            println("T : $t")
        end

        output[r] = sum(node_status)
    end

    return output

end

const N = 10^4
const threshold = fill(0.18, N)

@time data1 = diffusion_simulation(N, 50, threshold, 50, 100, "sync")
@time data2 = diffusion_simulation(N, 50, threshold, 50, 100, "async")
@time diffusion_simulation(N, 100, threshold, 100)


histogram(data1)
histogram(data2)

sum([i for i in data1 if i > 25])

histogram(degree(barabasi_albert(10000, 500)))
