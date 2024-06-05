# Author: Cameron Zachreson
# Institution: The University of Melbourne
# Simulation code acompanying the manuscript entitled: 
# "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
# Date released: Dec. 18, 2023


module Networks_RACF

using DataFrames
using StatsBase

using Main.Setup_RACF

abstract type Contact_T end 
abstract type Edge_T end 
abstract type E_list_T end 
abstract type E_list_temporal_T end 
abstract type N_list_T end 
abstract type N_list_temporal_multigraph_T end 


# contact structure for neighbour lists: 
mutable struct Contact <: Contact_T 
    id::Int64
    day::Int64
    weight::Float64
    room_id::Int64 #optional, not initialised by default 

    Contact(id, day, weight) = new(id, day, weight, 0)

end

mutable struct Edge <: Edge_T
    source_id::Int64
    target_id::Int64
    weight::Float64
end


mutable struct E_list <: E_list_T
    edges::Array{Edge_T, 1}
    E_list() = new(
        Array{Edge_T, 1}() # empty edge array. 
    )
end

#NOTE: this has a differrent structure to the temporal N_list
mutable struct E_list_temporal <: E_list_temporal_T
    day_to_E_list::Dict{Int64, E_list_T}
    
    E_list_temporal() = new(
        Dict{Int64, E_list_T}() #initialise an empty dict. day => E_list 
    )
end


mutable struct N_list <: N_list_T #neighbour list 
    id_to_contacts::Dict{Int64, Array{Contact_T, 1}}
    N_list() = new(
        Dict{Int64, Array{Contact_T, 1}}()
    )
end

mutable struct N_list_temporal_multigraph <: N_list_temporal_multigraph_T

    day_to_N_list::Dict{Int64, N_list_T} 
    # dict day => multigraph neighbour list 

    N_list_temporal_multigraph() = new(
        Dict{Int64, N_list_T}()
    )

end

#populates a dictionary source id => array[contact(target, day, weight)]
# see setup for dataframe read-in
function populate_neighbour_lists_from_DataFrame!(N_lists_out::Dict{Int64, Array{Contact_T, 1}}, N_lists_DF::DataFrame)

    for i = 1:size(N_lists_DF, 1)

        source_id = N_lists_DF.source_id[i]
    
        contacts_txt_i = split(N_lists_DF.out_edges[i], ';')
    
        n_contacts_i = size(contacts_txt_i, 1)
    
        for j = 1:n_contacts_i
    
            c_i_num = parse.(Int64, split(contacts_txt_i[j], ','))
            #format of edge is: [tartet id, day, weight]
            c_i = Contact(c_i_num[1], c_i_num[2], c_i_num[3])
    
            if haskey(N_lists_out, source_id)
                push!(N_lists_out[source_id], c_i)
            else
                N_lists_out[source_id] = [c_i]
            end
        end
    end
end


function sum_edge_weights_EList(elist::E_list_T)::Float64
    w_tot = 0.0
    for e in elist.edges
        w_tot += e.weight
    end
    return w_tot
end

function sum_edge_weights_contacts(contacts::Vector{Contact_T})::Float64
    w_tot = 0.0
    for c in contacts
        w_tot += c.weight 
    end
    return w_tot 
end


function compile_weights(elist::E_list_T)::Array{Float64}
    weights_out = []
    for e in elist.edges
        push!(weights_out, e.weight) #TODO check order. 
    end    

    return weights_out

end

function sample_E_list(elist::E_list_T, n, weight_arr, config::Setup_RACF.Config_T)::E_list_T

    e_list_out = E_list()

    edges = elist.edges

    edges_out = sample(config.rng_contacts, edges, Weights(weight_arr), n)

    e_list_out.edges = edges_out

    return e_list_out 

end

# note: this is time-consuming. scales as (n(d))^2
function compare_E_lists(elist_1::E_list_T, elist_2::E_list_T)::Bool

    lists_same = false
    n_dif = 0
    #n_dup = 0 # it's a multigraph, so we expect duplicates 
    n_same = 0
    n1 = size(elist_1.edges, 1)
    n2 = size(elist_2.edges, 1)
    if n1 != n2
        println("edge lists have different numbers of elements")
    else
        println("edge lists have the same number of elements")
    end


    # element-wise comparison. 
    # first, check that each element of elist_1 appears exactly once in elist_2
    # remember, default equality compares mutable structs by address 
    # see (e.g., https://stackoverflow.com/questions/70362843/how-to-create-equality-test-case-for-custom-structures-in-julia)
    for e1 in elist_1.edges
        n_e1 = 0
        for e2 in elist_2.edges
            # test equality by value 
            e1_e2 = edges_equal(e1, e2)
            if e1_e2
                n_e1 += 1
            end
        end

        if n_e1 == 0
            n_dif += 1
            println("edge ($e1) not found in edge list 2")
        else
            n_same += 1
        end
    end


    if (n1 == n2) && (n_same == n1)
        lists_same = true 
    end

    return lists_same 


end


# TODO: should generalise this for any mutable struct 
# so I can compare values in other cases. 
# can use the functions: 
#fields = fieldnames(typeof(s1))
#value = getfield(s1, fields[i])
# might have to be clever in how fields are matched comparing between structures. 
# e.g., could convert each struct to a dict, and compare key-value pairs. 
function edges_equal(e1::Edge_T, e2::Edge_T)::Bool
    
    equal_vals = false 

    if (e1.source_id === e2.source_id) 
        if (e1.target_id === e2.target_id)  
            if (e1.weight === e2.weight)
            equal_vals = true 
            end 
        end 
    end

    return equal_vals 

end


end