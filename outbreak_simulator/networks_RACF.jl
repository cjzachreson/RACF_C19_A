# contact structure for neighbour lists: 
mutable struct contact_type
    id::Int64
    day::Int64
    weight::Float64
end

mutable struct edge_type
    source_id::Int64
    target_id::Int64
    weight::Float64
end


mutable struct E_list
    edges::Array{edge_type}
    E_list() = new(
        [] # empty edge array. 
    )
end

#NOTE: this has a differrent structure to the temporal N_list
mutable struct E_list_temporal 
    day_to_E_list::Dict{Int64, E_list}
    
    E_list_temporal() = new(
        Dict{Int64, E_list}() #initialise an empty dict. day => E_list 
    )
end


mutable struct N_list_temporal

    id_to_contacts::Dict{Int64, Array{contact_type}}

    N_list_temporal() = new(
        Dict{Int64, Array{contact_type}}()
    )

end

#populates a dictionary source id => array[contact(target, day, weight)]
# see setup for dataframe read-in
function populate_neighbour_lists_from_DataFrame!(N_lists_out, N_lists_DF)

    for i = 1:size(N_lists_DF, 1)

        source_id = N_lists_DF.source_id[i]
    
        contacts_txt_i = split(N_lists_DF.out_edges[i], ';')
    
        n_contacts_i = size(contacts_txt_i, 1)
    
        for j = 1:n_contacts_i
    
            c_i_num = parse.(Int64, split(contacts_txt_i[j], ','))
            #format of edge is: [tartet id, day, weight]
            c_i = contact_type(c_i_num[1], c_i_num[2], c_i_num[3])
    
            if haskey(N_lists_out, source_id)
                push!(N_lists_out[source_id], c_i)
            else
                N_lists_out[source_id] = [c_i]
            end
        end
    end
end


function sum_edge_weights_EList(EL::E_list)::Float64
    w_tot = 0.0
    for e in EL.edges
        w_tot += e.weight
    end
    return w_tot
end

function sum_edge_weights_contacts(contacts::Vector{contact_type})::Float64
    w_tot = 0.0
    for c in contacts
        w_tot += c.weight 
    end
    return w_tot 
end


function compile_weights(EL::E_list)::Array{Float64}
    weights_out = []
    for e in EL.edges
        push!(weights_out, e.weight) #TODO check order. 
    end    

    return weights_out

end

function sample_E_list(EL::E_list, n, weight_arr)::E_list

    e_list_out = E_list()

    edges = EL.edges

    edges_out = sample(rng_infections, edges, Weights(weight_arr), n)

    e_list_out.edges = edges_out

    return e_list_out 

end