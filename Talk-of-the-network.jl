# *Note: This script is an effort to replicate the results from the paper "Talk
# of the Network: A Complex Systems Look at the Underlying Process of
# Word-of-Mouth", Goldenberg, Libai and Muller (2001). This is a self-didactic
# attempt.*
# For the writeup on the logic of this script please refer to the
# Jupyter notebook

using Distributions
srand(20130810)

mutable struct Node
    id::Int
    weak_ties::Vector{Int}
    strong_ties::Vector{Int}
    status::Bool
    activation_prob::Float64
end

function initialize_network(n_nodes::Int, n_strong_ties::Int, n_weak_ties::Int)

    # Initialize an empty network

    G = [Node(i, [], [], false, 0.0) for i in 1:n_nodes]
    node_ids = [node.id for node in G]

    # Wire the network according to the number of strong and weak ties
    # When wiring with random nodes, take care that the subject node and
    # already existing neighbors are not sampled again

    for node in G
        while length(node.weak_ties) < n_weak_ties
            rand_nbr = sample(node_ids[1:end .!= node.id])
            if !(rand_nbr in node.weak_ties || rand_nbr in node.strong_ties)
                push!(node.weak_ties, rand_nbr)
            end
        end
        while length(node.strong_ties) < n_strong_ties
            rand_nbr = sample(node_ids[1:end .!= node.id])
            if !(rand_nbr in node.weak_ties || rand_nbr in node.strong_ties)
                push!(node.strong_ties, rand_nbr)
            end
        end
    end

    return G
end

function reset_node_status!(G::Vector{Node})
    for node in G
        node.status = false
    end

    return nothing
end

function update_activation_prob!(G::Vector{Node}, alpha::Float64, beta_w::Float64, beta_s::Float64)

    for node in G

        n_active_weak_ties, n_active_strong_ties = 0, 0

        weak_ties = G[node.weak_ties]
        strong_ties = G[node.strong_ties]

        for weak_tie in weak_ties
            if weak_tie.status == true
                n_active_weak_ties += 1
            end
        end

        for strong_tie in strong_ties
            if strong_tie.status == true
                n_active_strong_ties += 1
            end
        end

        node.activation_prob = 1 - (1 - alpha) * (1 - beta_w)^n_active_weak_ties * (1 - beta_s)^n_active_strong_ties

    end

    return nothing
end

function update_status!(G::Vector{Node}, alpha::Float64, beta_w::Float64, beta_s::Float64)
    for node in G
        update_activation_prob!(G, alpha, beta_w, beta_s)

        if rand(Uniform()) < node.activation_prob
            node.status = true
        end
    end

    return nothing
end

function execute_simulation(parameter_space, n_nodes::Int, T::Int)

    # n_nodes dictates how big the network will be
    # T dictates the number of time steps for the simulation

    output = Matrix{Any}(length(parameter_space) * T, 7)
    colnames = ["s", "w", "alpha", "beta_w", "beta_s", "t", "num_engaged"]
    i = 1 # Index the output

    # Rewiring the network each time is expensive. We can cut down repeats of the same rewiring process
    # by building the network only when the parameters used to build the network have changed.

    old_s, old_w = parameter_space[1][1:2]
    G = initialize_network(n_nodes, old_s, old_w)

    for (s, w, alpha, beta_w, beta_s) in parameter_space[1:end]

        # Rewire the network only if the network creation parameters have changed
        # Once initialize_network is fast, this portion of the code can be removed

        if !(old_s == s && old_w == w)
            G = initialize_network(n_nodes, s, w)
        end
        reset_node_status!(G)

        println("Beginning simulation on setting $((s, w, alpha, beta_w, beta_s)) at : ", Dates.format(now(), "HH:MM"))
        for t in 1:T
            update_status!(G, alpha, beta_w, beta_s)
            num_engaged = sum([node.status for node in G])
            output[i, :] = [s, w, alpha, beta_w, beta_s, t, num_engaged]
            i += 1
        end

        old_s, old_w = s, w
    end

    return [colnames; output]
end

parameter_space = [(s, w, alpha, beta_w, beta_s) for s in floor.(Int, linspace(5, 29, 7)),
                                                     w in floor.(Int, linspace(5, 29, 7)),
                                                     alpha in linspace(0.0005, 0.01, 7),
                                                     beta_w in linspace(0.005, 0.015, 7),
                                                     beta_s in linspace(0.01, 0.07, 7)]

@time results = execute_simulation(parameter_space, 3000, 30)


writecsv("test.csv", 1:1000)
