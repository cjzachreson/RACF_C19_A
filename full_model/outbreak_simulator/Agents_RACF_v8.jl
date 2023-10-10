#agent structure type definiition: 
# see setup for abstract type definitions

module Agents_RACF

using DataFrames 
using Distributions

using Main.Setup_RACF
using Main.Networks_RACF
using Main.Diseases_RACF

abstract type Agent_T end
abstract type Agents_T end 

mutable struct Worker_Agent <: Agent_T
    id::Int64 #unique agent ID
    is_medical::Bool # Boolean flag, true means this is a medical worker
    roster::Array{Int64, 1} # Vector length is number of days over which schedule repeats 0s  1s 
    rooms::Dict{Int64, Array{Int64}} #maps days i.e., 1-7 to array of room ids worked on each day
    n_rooms::Int64 # number of room-days serviced by this worker over a roster period
    
    ## ubiquitous parameters 
    contacts::Dict{Int64, Array{Networks_RACF.Contact_T}} # maps day -> [target_id, weight, day]

    infections::Dict{String, Diseases_RACF.Infection_T} # pathogen name  => infection
    n_infections::Int64 #total number of times infected (all time)
    immunity::Diseases_RACF.Immunity_Profile_T # agent's immunity profile 
    t_detected::Float64 # time of most recent positive test
    t_removed::Float64 # positive if agent is removed from the facility (reversible)

    # testing schedule, same structure as roster
    test_schedule::Array{Int64, 1} # Vector length is number of days over which schedule repeats 0s  1s 

    #default initialiser: 
    Worker_Agent(id) = new(id, #id::Int64
    false, #is_medical::Bool
    Array{Int64, 1}(), #roster::Array{Int64}
    Dict{Int64, Array{Int64}}(), # rooms, instantiated as Dict, day => rooms. 
    0, #n_rooms::Int64
    ## ubiquitous agent parameteers 
    Dict{Int64, Array{Networks_RACF.Contact_T}}(), #contacts::Dict{Int64, Array{contact_type}} #Dict allows indexing by day
    Dict{String, Diseases_RACF.Infection_T}(),# empty infection dict 
    0, #number of infections so far 
    Diseases_RACF.Immunity_Profile(), #initialise an empty immunity profile. 
    -1.0, # initialise test time to negative (no tests yet.)
    -1.0, #t_removed
    Array{Int64, 1}() # test_schedule 
    ) 
end


mutable struct Resident_Agent <: Agent_T
    id::Int64 #unique agent ID
    chc_needs::Bool # Boolean flag, true means this resident has high medical needs
    beh_needs::Bool # Boolean flag, true means this resident has high behavioural needs
    adl_needs::Bool # Boolean flag, true means this resident has high daily living needs
    room::Int64 # resident's room id
    
    contacts::Dict{Int64, Array{Networks_RACF.Contact_T}} # maps day -> [target_id, day, weight]

    infections::Dict{String, Diseases_RACF.Infection_T} # pathogen name  => infection
    n_infections::Int64 #total number of times infected (all time)
    immunity::Diseases_RACF.Immunity_Profile_T # agent's immunity profile 
    t_detected::Float64 # time of most recent positive test
    t_removed::Float64 # positive if agent is removed from the facility (reversible)


    #default initialiser: 
    Resident_Agent(id) = new(id, #id::Int64
    false, #chc_needs::Bool # Boolean flag, true means this resident has high medical needs
    false, #beh_needs::Bool # Boolean flag, true means this resident has high behavioural needs
    false, #adl_needs::Bool # Boolean flag, true means this resident has high daily living needs
    0, #room::Int64 # resident's room id
    Dict{Int64, Array{Networks_RACF.Contact_T}}(), #contacts::Dict{Int64, Array{contact_type}} # maps day -> [target_id, day, weight]
    Dict{String, Diseases_RACF.Infection_T}(),# empty infection dict 
    0, #number of infections so far 
    Diseases_RACF.Immunity_Profile(), #initialise an empty immunity profile. 
    -1.0, # initialise test time to negative (no tests yet.)
    -1.0 #t_removed
    ) 
end


mutable struct Agents <: Agents_T

    workers_G::Dict{Int64, Agent_T}
    workers_M::Dict{Int64, Agent_T}
    residents::Dict{Int64, Agent_T}
    All::Dict{Int64, Agent_T}

    Agents() = new(
        Dict{Int64, Agent_T}(),
        Dict{Int64, Agent_T}(),
        Dict{Int64, Agent_T}(),
        Dict{Int64, Agent_T}()
    )

end



# functions for property assignment from dataframe:

# adds immunity status


#adds immunity status from distribution, rather than from file: 
function initialise_immunity_status_from_pdist!(im::Diseases_RACF.Immunity_Profile_T, 
                                                pdist::Distribution, 
                                                config::Setup_RACF.Config_T)

    #input is immunity status of agent i, and 
    # text dataframe of immunity status values from population
    # generator output (only agent i's line of the DF)

    # NOTE: still treating neuts as static for now, 


    im.NAT_peak = 0.0

    im.dt_peak = 0.0

    im.NAT_decay_rate = 0.0

    im.NAT_t = 0.0
 

    if config.immunity
        # sample from population-specific neut distribution 
        ln_neuts = rand(config.rng_immunity, pdist)

        im.protection_Infection["Default"] = Diseases_RACF.Efficacy_infection_Default(ln_neuts)

        im.protection_Symptoms["Default"] = Diseases_RACF.Efficacy_Symptoms_Default(ln_neuts)

        im.protection_Transmission["Default"]  = Diseases_RACF.Efficacy_OT_Default(ln_neuts)

        im.protection_Death["Default"] = Diseases_RACF.Efficacy_Death_Default(ln_neuts)
    else

        im.protection_Infection["Default"] = 0.0

        im.protection_Symptoms["Default"] = 0.0

        im.protection_Transmission["Default"]  = 0.0

        im.protection_Death["Default"] = 0.0
    end


    #Check that data types are parsed correctly. 

end

function add_neighbours_from_N_list_v2!(a::Agent_T, 
                                        d::Int64, 
                                        id_to_contacts::Dict{Int64, Array{Networks_RACF.Contact_T, 1}},
                                        id_i::Int64)
    #parse contacts and populate: 
    if haskey(id_to_contacts, id_i) 
        #normally not necessary to check the id_i key, 
        #as N_lists will compile all days of the roster - if the worker is not 
        # in the edge list, something is wrong

        contacts_i = id_to_contacts[id_i]
        
        a.contacts[d] = contacts_i 
    end
end

function assign_contacts!(agents::Agents_T, 
                          N_lists_TMG::Networks_RACF.N_list_temporal_multigraph_T)

    for (d, N_list) in N_lists_TMG.day_to_N_list

        id_to_contacts = N_list.id_to_contacts
        for (id_i, a) in agents.All
            add_neighbours_from_N_list_v2!(a, d, id_to_contacts, id_i)
        end

    end
end

function assign_todays_contacts!(agents::Agents_T, 
                                 N_lists_TMG::Networks_RACF.N_list_temporal_multigraph_T, 
                                 day_of_week::Int64) 

    N_list = N_lists_TMG.day_to_N_list[day_of_week]
    id_to_contacts = N_list.id_to_contacts
    for (id_i, a) in agents.All # in updated version this will operate only on a subset
        add_neighbours_from_N_list_v2!(a, day_of_week, id_to_contacts, id_i)
    end

end
# 2022 09 14, removed initialisation of N_lists from text file.
#populates a dictionary of worker agents [id -> agent object]
function populate_workers_from_DataFrame!(workers_out::Dict{Int64, Agent_T}, 
                                          workers_DF::DataFrame,
                                          config::Setup_RACF.Config_T)#, N_lists)  
    
    n_workers = size(workers_DF, 1)
    
    for i = 1:n_workers
    
        id_i =  workers_DF.id[i]
    
        workers_out[id_i] = Worker_Agent(id_i)
    
        workers_out[id_i].is_medical = workers_DF.medical[i]
        # parse room string (parse implementation)
        workers_out[id_i].roster = parse.(Int64, split(workers_DF.roster[i], ','))
    
        #parse rooms and populate: 
        rooms_txt_i = split(workers_DF.rooms[i], ';') #makes a vector with an element for each room
        
        n_rooms_i = size(rooms_txt_i, 1)
        
        for r = 1:n_rooms_i
            r_num = parse.(Int64, split(rooms_txt_i[r], ','))
            d_r = r_num[1]
            rm_id_r = r_num[2]
            if haskey(workers_out[id_i].rooms, d_r) #useful - after construction, no key means no rooms
                push!(workers_out[id_i].rooms[d_r], rm_id_r)
            else
                workers_out[id_i].rooms[d_r] = [rm_id_r]
            end
        end
    
        #add_neighbours_from_N_list!(workers_out[id_i], N_lists, id_i)
    
        initialise_immunity_status!(workers_out[id_i].immunity, 
                                    workers_DF[workers_DF.id .== id_i, :], 
                                    config)
    
    end

end

# 2022 09 14, removed initialisation of N_lists from text file.
#populates a dictionary of resident agents [id -> agent object]
function populate_residents_from_DataFrame!(residents_out::Dict{Int64, Agent_T}, 
                                            residents_DF::DataFrame,
                                            config::Setup_RACF.Config_T)#, N_lists)
    
    n_residents = size(residents_DF, 1)
    
    for i = 1:n_residents
        
        #identifer 
        id_i =  residents_DF.id[i]
    
        # initialise
        residents_out[id_i] = Resident_Agent(id_i)
    
        # needs
        residents_out[id_i].chc_needs = residents_DF.chc_needs[i]
        residents_out[id_i].beh_needs = residents_DF.beh_needs[i]
        residents_out[id_i].adl_needs = residents_DF.adl_needs[i]
        
        # room
        residents_out[id_i].room = residents_DF.room_id[i]

        #parse contacts and populate: 
        #add_neighbours_from_N_list!(residents_out[id_i], N_lists, id_i)

        #add immunity status: 
        initialise_immunity_status!(residents_out[id_i].immunity, 
                                    residents_DF[residents_DF.id .== id_i, :], 
                                    config)

    end
end

#2022 09 25, new version of populate_residents that takes immunity status from a distribution 
# instead of taking it from the input file. 
function populate_workers_from_DataFrame_imDist!(workers_out::Dict{Int64, Agent_T}, 
                                                 workers_DF::DataFrame, 
                                                 neut_dist::Distribution,
                                                 config::Setup_RACF.Config_T)#, N_lists)  
    
    n_workers = size(workers_DF, 1)
    
    for i = 1:n_workers
    
        id_i =  workers_DF.id[i]
    
        workers_out[id_i] = Worker_Agent(id_i)
    
        workers_out[id_i].is_medical = workers_DF.medical[i]
        # parse room string (parse implementation)
        workers_out[id_i].roster = parse.(Int64, split(workers_DF.roster[i], ','))
    
        #parse rooms and populate: 
        rooms_txt_i = split(workers_DF.rooms[i], ';') #makes a vector with an element for each room
        
        n_rooms_i = size(rooms_txt_i, 1)
        
        for r = 1:n_rooms_i
            r_num = parse.(Int64, split(rooms_txt_i[r], ','))
            d_r = r_num[1]
            rm_id_r = r_num[2]
            if haskey(workers_out[id_i].rooms, d_r) #useful - after construction, no key means no rooms
                push!(workers_out[id_i].rooms[d_r], rm_id_r)
            else
                workers_out[id_i].rooms[d_r] = [rm_id_r]
            end
        end
    
        #add_neighbours_from_N_list!(workers_out[id_i], N_lists, id_i)
    
        initialise_immunity_status_from_pdist!(workers_out[id_i].immunity, neut_dist, config)
    
    end

end

function populate_residents_from_DataFrame_imDist!(residents_out::Dict{Int64, Agent_T}, 
                                                   residents_DF::DataFrame, 
                                                   neut_dist::Distribution,
                                                   config::Setup_RACF.Config_T)#, N_lists)
    
    n_residents = size(residents_DF, 1)
    
    for i = 1:n_residents
        
        #identifer 
        id_i =  residents_DF.id[i]
    
        # initialise
        residents_out[id_i] = Resident_Agent(id_i)
    
        # needs
        residents_out[id_i].chc_needs = residents_DF.chc_needs[i]
        residents_out[id_i].beh_needs = residents_DF.beh_needs[i]
        residents_out[id_i].adl_needs = residents_DF.adl_needs[i]
        
        # room
        residents_out[id_i].room = residents_DF.room_id[i]

        #parse contacts and populate: 
        #add_neighbours_from_N_list!(residents_out[id_i], N_lists, id_i)

        #add immunity status: 
        initialise_immunity_status_from_pdist!(residents_out[id_i].immunity, neut_dist, config)

    end
end

#utilities: 
function is_worker(a::Agent_T)::Bool
    return (typeof(a) == Worker_Agent)
end

function is_resident(a::Agent_T)::Bool
    return (typeof(a) == Resident_Agent)
end

function is_high_needs(a::Agent_T)::Bool
    if is_resident(a)
    return (a.adl_needs || a.beh_needs || a.chc_needs)
    else 
        return false 
    end
end

function count_high_needs_residents(agents::Agents_T)::Int64
    n = 0
    for (id, a) in agents.residents
        if is_high_needs(a)
            n += 1
        end
    end
    return n 
end 

function count_reg_needs_residents(agents::Agents_T)::Int64
    n = 0
    for (id, a) in agents.residents
        if !(is_high_needs(a))
            n += 1
        end
    end
    return n 
end

# model-specific pairwise contact weighting (i.e., high needs etc.)
# NOTE: this function modifies the weight property of each contact
function pairwise_weights!(A::Agents_T, config::Setup_RACF.Config_T)
    for (id, source) in A.All
        c_set = source.contacts # Dict::{day::int64, contacts::Array{contact_type}}
        # iterate through days
        for (d, c) in c_set
            #iterate through the array of contacts on day d
            for c_i in c 
                target = A.All[c_i.id]
                # check type of source and target
                # TODO: functions like this should be generalised based on some
                # user-specified mapping. [source type, target type] => [weight factor] 
                w = 1.0 #default pairwise factor 
                if (is_resident(source) && is_worker(target))
                    if is_high_needs(source)
                        w = config.w_high_needs # GLOBAL
                    else 
                        w = config.w_reg_needs 
                    end
                elseif (is_worker(source) && is_resident(target))
                    if is_high_needs(target)
                        w = config.w_high_needs # GLOBAL
                    else 
                        w = config.w_reg_needs 
                    end
                elseif (is_worker(source) && is_worker(target))
                    w = config.w_worker_worker
                elseif (is_resident(source) && is_resident(target))
                    if source.room == target.room
                        w = config.w_res_same_room
                    else
                        w = config.w_res_diff_room
                    end
                end

                #c_i.weight *= w # using multigraph now, and updating weight each day
                c_i.weight = w 
            
            end
        end
    end
end

function pairwise_weights_d!(A::Agents_T, day_of_week::Int64, config::Setup_RACF.Config_T)
    for (id, source) in A.All
        c_set = source.contacts # Dict::{day::int64, contacts::Array{contact_type}}
        # iterate through days
        #for (d, c) in c_set
        if haskey(c_set, day_of_week)
            c = c_set[day_of_week]
            #iterate through the array of contacts on day d
            for c_i in c 
                target = A.All[c_i.id]
                # check type of source and target
                # TODO: functions like this should be generalised based on some
                # user-specified mapping. [source type, target type] => [weight factor] 
                w = 1.0 #default pairwise factor 
                if (is_resident(source) && is_worker(target))
                    if is_high_needs(source)
                        w = config.w_high_needs # GLOBAL
                    else 
                        w = config.w_reg_needs 
                    end
                elseif (is_worker(source) && is_resident(target))
                    if is_high_needs(target)
                        w = config.w_high_needs # GLOBAL
                    else 
                        w = config.w_reg_needs 
                    end
                elseif (is_worker(source) && is_worker(target))
                    w = config.w_worker_worker
                elseif (is_resident(source) && is_resident(target))
                    if source.room == target.room
                        w = config.w_res_same_room
                    else
                        w = config.w_res_diff_room
                    end
                end

                c_i.weight = w 
            
            end
        end
    end
end


function add_source_edges_to_E_list!(elist::Networks_RACF.E_list_T, 
                                     a::Agent_T, 
                                     contacts::Array{Networks_RACF.Contact_T})
    source_id = a.id
    
    for c in contacts 
        target_id = c.id 
        weight = c.weight 
        e = Networks_RACF.Edge(source_id, target_id, weight) 

        push!(elist.edges, e)
    end

end


function fill_E_lists_all!(elist_t::Networks_RACF.E_list_temporal_T, 
                           A::Agents_T)

    # iterate through agents and contacts
    for (id, a) in A.All
        for (d, c) in a.contacts
            #add edges to lists 
            if haskey(elist_t.day_to_E_list, d)

            add_source_edges_to_E_list!(elist_t.day_to_E_list[d], a, c)
            
            else # initialise for day d  
                elist_t.day_to_E_list[d] = Networks_RACF.E_list()
                add_source_edges_to_E_list!(elist_t.day_to_E_list[d], a, c)
            end
        end
    end
end


function compute_total_weight_d(a::Agent_T, d::Int64)::Float64

    w = 0.0
    for c in a.contacts[d]
        w += c.weight
    end

    return w 

end


function compute_total_FTE(agents::Agents_T)::Float64

    total_FTE = 0
    for (id, a) in agents.All
        if is_worker(a)
            total_FTE += (sum(a.roster) / 5)
        end
    end

    return total_FTE
end

function compute_FTE_deficit(agents::Agents_T, removed_workers::Dict{Int64, Float64})::Float64

    FTE_deficit = 0 

    for (id, t) in removed_workers
        FTE_deficit += (sum(agents.All[id].roster) / 5)
    end

    return FTE_deficit 

end


function count_worker_infections(agents::Agents_T, all_transmissions::DataFrame)::Int64
    count = 0
    for id in all_transmissions.target
        a = agents.All[id]
        if is_worker(a)
            count += 1
        end
    end
    return count
end

function count_resident_infections(agents::Agents_T, all_transmissions::DataFrame)::Int64
    count = 0
    for id in all_transmissions.target
        a = agents.All[id]
        if is_resident(a)
            count += 1
        end
    end
    return count

end

function count_worker_detections(agents::Agents_T, all_detections::DataFrame)::Int64
    count = 0
    for id in all_detections.id
        a = agents.All[id]
        if is_worker(a)
            count += 1
        end
    end
    return count

end

function count_resident_detections(agents::Agents_T, all_detections::DataFrame)::Int64
    count = 0
    for id in all_detections.id
        a = agents.All[id]
        if is_resident(a)
            count += 1
        end
    end
    return count
end




end