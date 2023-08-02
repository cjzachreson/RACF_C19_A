#agent structure type definiition: 
# see setup for abstract type definitions

mutable struct worker_agent_type <: Agent_T
    id::Int64 #unique agent ID
    is_medical::Bool # Boolean flag, true means this is a medical worker
    roster::Array{Int64} # Vector length is number of days over which schedule repeats 0s  1s 
    rooms::Dict{Int64, Array{Int64}} #maps days i.e., 1-7 to array of rooms worked on each day
    n_rooms::Int64 # number of room-days serviced by this worker over a roster period
    
    ## ubiquitous parameters 
    contacts::Dict{Int64, Array{contact_type}} # maps day -> [target_id, weight, day]

    infections::Dict{String, Infection_T} # pathogen name  => infection
    n_infections::Int64 #total number of times infected (all time)
    immunity::immunity_profile # agent's immunity profile 
    t_detected::Float64 # time of most recent positive test
    t_removed::Float64 # positive if agent is removed from the facility (reversible)

    #default initialiser: 
    worker_agent_type(id) = new(id, #id::Int64
    false, #is_medical::Bool
    [], #roster::Array{Int64}
    Dict{Int64, Array{Int64}}(), # rooms, instantiated as Dict, day => rooms. 
    0, #n_rooms::Int64
    ## ubiquitous agent parameteers 
    Dict{Int64, Array{contact_type}}(), #contacts::Dict{Int64, Array{contact_type}} #Dict allows indexing by day
    Dict{String, Infection_T}(),# empty infection dict 
    0, #number of infections so far 
    immunity_profile(), #initialise an empty immunity profile. 
    -1.0, # initialise test time to negative (no tests yet.)
    -1.0 #t_removed 
    ) 
end


mutable struct resident_agent_type <: Agent_T
    id::Int64 #unique agent ID
    chc_needs::Bool # Boolean flag, true means this resident has high medical needs
    beh_needs::Bool # Boolean flag, true means this resident has high behavioural needs
    adl_needs::Bool # Boolean flag, true means this resident has high daily living needs
    room::Int64 # resident's room id
    
    contacts::Dict{Int64, Array{contact_type}} # maps day -> [target_id, day, weight]

    infections::Dict{String, Infection_T} # pathogen name  => infection
    n_infections::Int64 #total number of times infected (all time)
    immunity::immunity_profile # agent's immunity profile 
    t_detected::Float64 # time of most recent positive test
    t_removed::Float64 # positive if agent is removed from the facility (reversible)

    #default initialiser: 
    resident_agent_type(id) = new(id, #id::Int64
    false, #chc_needs::Bool # Boolean flag, true means this resident has high medical needs
    false, #beh_needs::Bool # Boolean flag, true means this resident has high behavioural needs
    false, #adl_needs::Bool # Boolean flag, true means this resident has high daily living needs
    0, #room::Int64 # resident's room id
    Dict{Int64, Array{contact_type}}(), #contacts::Dict{Int64, Array{contact_type}} # maps day -> [target_id, day, weight]
    Dict{String, Infection_T}(),# empty infection dict 
    0, #number of infections so far 
    immunity_profile(), #initialise an empty immunity profile. 
    -1.0, # initialise test time to negative (no tests yet.)
    -1.0 #t_removed 
    ) 
end


mutable struct Agents

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
function initialise_immunity_status!(im::immunity_profile, agent_DF)

    #input is immunity status of agent i, and 
    # text dataframe of immunity status values from population
    # generator output (only agent i's line of the DF)

    # TODO: convert numeric week to datetime and add vaccine and 
    # infection histories. (this will require a map from numeric week to datetime.)
    # vaccination history Dict{DateTime, String}
    # infection history Dict{DateTime, Disease_T}

    #NAT_peak::Float64
    im.NAT_peak = 0.0
    #dt_peak::Float64
    im.dt_peak = 0.0
    #NAT_decay_rate::Float64
    im.NAT_decay_rate = 0.0
    #NAT_t::Float64
    im.NAT_t = 0.0
    # protection values 

    if IMMUNITY
        #for now, these are all for Omicron (SARS-CoV-2)
        #protection_Infection::Dict{String, Float64}
        im.protection_Infection["Omicron"] = agent_DF.eff_acqui[1]
        #protection_Symptoms::Dict{String, Float64}
        im.protection_Symptoms["Omicron"] = agent_DF.eff_symp[1]
        #protection_Transmission::Dict{String, Float64}
        im.protection_Transmission["Omicron"]  = agent_DF.eff_trans[1]
        #protection_Death::Dict{String, Float64} 
        im.protection_Death["Omicron"] = agent_DF.eff_death[1]
    else
        #for now, these are all for Omicron (SARS-CoV-2)
        #protection_Infection::Dict{String, Float64}
        im.protection_Infection["Omicron"] = 0.0
        #protection_Symptoms::Dict{String, Float64}
        im.protection_Symptoms["Omicron"] = 0.0
        #protection_Transmission::Dict{String, Float64}
        im.protection_Transmission["Omicron"]  = 0.0
        #protection_Death::Dict{String, Float64} 
        im.protection_Death["Omicron"] = 0.0
    end


    #Check that data types are parsed correctly. 

end

# adds contact structure (bi-directional, treated here as out-edges) 
function add_neighbours_from_N_list(a::Agent_T, N_lists, id_i)
        #parse contacts and populate: 
        if haskey(N_lists, id_i) 
            #normally not necessary to check the id_i key, 
            #as N_lists will compile all days of the roster - if the worker is not 
            # in the edge list, something is wrong
    
            contacts_i = N_lists[id_i]
            
            n_contacts_i = length(contacts_i)
    
            for j = 1:n_contacts_i
    
                day = contacts_i[j].day
                #target_id = contacts_i[j].id
                #weight = contacts_i[j].weight
    
                if haskey(a.contacts, day)
                push!(a.contacts[day], contacts_i[j])
                else
                    #appears to pass a copy rather than a reference
                    # this means N_lists can be modified without modifying these. 
                    a.contacts[day] = [contacts_i[j]] 
                end 
            end
        end
end

#populates a dictionary of worker agents [id -> agent object]
function populate_workers_from_DataFrame!(workers_out::Dict{Int64, Agent_T}, workers_DF, N_lists)
    
    n_workers = size(workers_DF, 1)
    #workers_out = Dict{Int64, worker_agent_type}()
    
    for i = 1:n_workers
    
        id_i =  workers_DF.id[i]
    
        workers_out[id_i] = worker_agent_type(id_i)
    
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
    
        add_neighbours_from_N_list(workers_out[id_i], N_lists, id_i)
    
        initialise_immunity_status!(workers_out[id_i].immunity, workers_DF[workers_DF.id .== id_i, :])
    
    end

end

#populates a dictionary of resident agents [id -> agent object]
function populate_residents_from_DataFrame!(residents_out::Dict{Int64, Agent_T}, residents_DF, N_lists)
    
    n_residents = size(residents_DF, 1)
    
    for i = 1:n_residents
        
        #identifer 
        id_i =  residents_DF.id[i]
    
        # initialise
        residents_out[id_i] = resident_agent_type(id_i)
    
        # needs
        residents_out[id_i].chc_needs = residents_DF.chc_needs[i]
        residents_out[id_i].beh_needs = residents_DF.beh_needs[i]
        residents_out[id_i].adl_needs = residents_DF.adl_needs[i]
        
        # room
        residents_out[id_i].room = residents_DF.room_id[i]

        #parse contacts and populate: 
        add_neighbours_from_N_list(residents_out[id_i], N_lists, id_i)

        #add immunity status: 
        initialise_immunity_status!(residents_out[id_i].immunity, residents_DF[residents_DF.id .== id_i, :])

    end
end


#utilities: 
function is_worker(a::Agent_T)::Bool
    return (typeof(a) == worker_agent_type)
end

function is_resident(a::Agent_T)::Bool
    return (typeof(a) == resident_agent_type)
end

function is_high_needs(a::resident_agent_type)::Bool
    return (a.adl_needs || a.beh_needs || a.chc_needs)
end

# model-specific pairwise contact weighting (i.e., high needs etc.)
# NOTE: this function modifies the weight property of each contact
function pairwise_weights!(A::Agents)
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
                        w = w_high_needs # GLOBAL
                    else 
                        w = w_reg_needs 
                    end
                elseif (is_worker(source) && is_resident(target))
                    if is_high_needs(target)
                        w = w_high_needs # GLOBAL
                    else 
                        w = w_reg_needs 
                    end
                elseif (is_worker(source) && is_worker(target))
                    w = w_worker_worker
                elseif (is_resident(source) && is_resident(target))
                    if source.room == target.room
                        w = w_res_same_room
                    else
                        w = w_res_diff_room
                    end
                end

                c_i.weight *= w 
            
            
            
            
            end
        
        
        
        
        end
    end
end


function add_source_edges_to_E_list!(EL::E_list, a::Agent_T, contacts::Array{contact_type})
    source_id = a.id
    
    for c in contacts 
        target_id = c.id 
        weight = c.weight 
        e = edge_type(source_id, target_id, weight) 

        push!(EL.edges, e)
    end

end


function fill_E_lists_all!(EL_t::E_list_temporal, A::Agents)

    # iterate through agents and contacts
    for (id, a) in A.All
        for (d, c) in a.contacts
            #add edges to lists 
            if haskey(EL_t.day_to_E_list, d)

            add_source_edges_to_E_list!(EL_t.day_to_E_list[d], a, c)
            
            else # initialise for day d  
                EL_t.day_to_E_list[d] = E_list()
                add_source_edges_to_E_list!(EL_t.day_to_E_list[d], a, c)
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







# recovery


