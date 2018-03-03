using Distributions, LightGraphs
using DataFrames

mutable struct ThresholdModel
    G::LightGraphs.SimpleGraphs.SimpleGraph{Int}
    node_status::BitVector
    threshold::Vector{T} where T <: Real
    threshold_type::String
    t::Int
    n_realizations::Int
end

# Watts models
function fraction_engaged_nbrs(tm::ThresholdModel, node::Int)
    """
    Computes the fraction of neighbors engaged within the neighborhood
    of a given node. It uses the node status to check the engagement status of
    each of the nodes neighbors.

    Useful for models with fractional thresholds
    """
    n_engaged_nbrs = sum([tm.node_status[nbr] for nbr in neighbors(tm.G, node)])

    return n_engaged_nbrs/length(neighbors(tm.G, node))
end

# Centola-Macy models
function num_engaged_nbrs(tm::ThresholdModel, node::Int)
    """
    Computes the number of neighbors engaged within the neighborhood
    of a given node. It uses the node status to check the engagement status of
    each of the nodes neighbors.

    Useful for models with integer thresholds
    """
    return sum([tm.node_status[nbr] for nbr in neighbors(tm.G, node)])
end

function update_node_status!(tm::ThresholdModel)
    """
    This function executes the random asynchronous updates of the entire network
    at each time step. In this conceptualization, each time step comprises mini
    time steps during which a randomly shuffled node list updates.

    The node status vector is updated in place.
    """

    if tm.threshold_type == "frac"
        for node in shuffle(vertices(tm.G))
            if tm.node_status[node] == false
                if fraction_engaged_nbrs(tm, node) > tm.threshold[node]
                    tm.node_status[node] = true
                end
            end
        end
    elseif tm.threshold_type == "int"
        for node in shuffle(vertices(tm.G))
            if tm.node_status[node] == false
                if num_engaged_nbrs(tm, node) > tm.threshold[node]
                    tm.node_status[node] = true
                end
            end
        end
    end

    return nothing
end
