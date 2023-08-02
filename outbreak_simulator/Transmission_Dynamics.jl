

# infectious disease dynamics (functions)

# includes functions for: 

# infection: 

function infect_agent!(a::Agent_T, pathogen::Disease_T, time_of_exposure::Float64, infected_agents::Dict{Int64, Float64})

    t_inc = rand(rng_infections, pathogen.inc_dist)
    q_inc = cdf(pathogen.inc_dist, t_inc)
    t_rec = rand(rng_infections, pathogen.rec_dist)
    beta_max = rand(rng_infections, pathogen.b_dist)

    a.n_infections += 1
    new_infection = infection(pathogen, 0.0, t_inc, q_inc, t_rec, beta_max) # time since infection is 0
    
    # reinitialise symptom expression based on immunity status: 
    p_symp_i = 1.0 - pathogen.p_asymp
    if haskey(a.immunity.protection_Symptoms, new_infection.pathogen_name)
        eff_symp = 1.0 - a.immunity.protection_Symptoms[new_infection.pathogen_name]
        p_symp_i *= eff_symp
    end

    if rand(rng_infections) < p_symp_i
        new_infection.symptomatic = true
    else
        new_infection.symptomatic = false
    end
    

    a.infections[pathogen.name] = new_infection

    infected_agents[a.id] = time_of_exposure #updating dictionary of infected agents 
end

function is_infected(a::Agent_T, p_name::String)::Bool
    if haskey(a.infections, p_name)
        return true
    else
        return false
    end
end

# pairwise multistrain transmission: 
# TODO: think about pairwise vs. target-integrated transmission... 
function transmit_infection!(source::Agent_T, target::Agent_T, time_of_exposure::Float64, infected_agents::Dict{Int64, Float64})::Bool

    infection_occurred = false

    for (key, infection) in source.infections 

        #compute a probability of transmission 
        pathogen_name = infection.pathogen_name
        # break if target is already infected with the same pathogen: 

        if !is_infected(target, pathogen_name)

            s_im = source.immunity
            t_im = target.immunity 

            if t_im.protection_Infection[pathogen_name] < 1.0 # target not totally immune

                foi_source = infection.beta_t
                immunity_fac_source = 0.0
                immunity_fac_target = 0.0 

                if haskey(s_im.protection_Transmission, pathogen_name)
                    immunity_fac_source = 1.0 - s_im.protection_Transmission[pathogen_name]
                end
                if haskey(t_im.protection_Infection, pathogen_name)
                    immunity_fac_target = 1.0 - t_im.protection_Infection[pathogen_name]
                end

                foi_source *= immunity_fac_source
                foi_source *= immunity_fac_target 

                # TODO: include other factors associated with transmission. 
                p_trans = 1.0 - exp((-1.0 * foi_source * dt))

                #println("agent $(source.id) trying to infect agent $(target.id) with $pathogen_name with probability $p_trans, b_max = $(source.infections[pathogen_name].beta_max)")

                #if p_trans == 0.0 && infection.t_infected > infection.t_latent 
                 #   println("* PROBLEM with FOI compuation  *")
                #end
                
                if rand(rng_infections) < p_trans
                    infect_agent!(target, infection.pathogen, time_of_exposure, infected_agents)

                    #if target.id == 119
                     #   println("**agent $(source.id) infected agent $(target.id) with $pathogen_name with probability $p_trans")
                    #end

                    infection_occurred =  true
                end
            end
        end
    end 

    return infection_occurred
end

function transmit_infection_AOB!(source::Agent_T, target::Agent_T, time_of_exposure::Float64, infected_agents::Dict{Int64, Float64})::Bool

    infection_occurred = false

    for (key, infection) in source.infections 

        #compute a probability of transmission 
        pathogen_name = infection.pathogen_name
        # break if target is already infected with the same pathogen: 

        if !is_infected(target, pathogen_name)

            s_im = source.immunity
            t_im = target.immunity 

            st_IC = 0.0 
            # compute reduction in FoI based on Infection Control measurs used during active outbreaks: 
            if ( is_worker(target) && is_worker(source) )
                st_IC = eff_IC_worker_worker #global 
            elseif ( (is_worker(target) && is_resident(source)) || (is_worker(source) && is_resident(target) ))
                st_IC = eff_IC_worker_resident # global 
            elseif ( is_resident(target) && is_resident(source) )
                st_IC = eff_IC_resident_resident # global 
            end 

            IC_fac = 1.0 - st_IC

            if t_im.protection_Infection[pathogen_name] < 1.0 # target not totally immune

                foi_source = infection.beta_t
                immunity_fac_source = 0.0
                immunity_fac_target = 0.0 

                if haskey(s_im.protection_Transmission, pathogen_name)
                    immunity_fac_source = 1.0 - s_im.protection_Transmission[pathogen_name]
                end
                if haskey(t_im.protection_Infection, pathogen_name)
                    immunity_fac_target = 1.0 - t_im.protection_Infection[pathogen_name]
                end

                foi_source *= immunity_fac_source 
                foi_source *= immunity_fac_target
                foi_source *= IC_fac # infection control  

                # TODO: include other factors associated with transmission. 
                p_trans = 1.0 - exp((-1.0 * foi_source * dt))

                #println("agent $(source.id) trying to infect agent $(target.id) with $pathogen_name with probability $p_trans")

                
                #if p_trans == 0.0 && infection.t_infected > infection.t_latent 
                 #   println("* PROBLEM with FOI compuation  *")
                #end
                
                if rand(rng_infections) < p_trans
                    infect_agent!(target, infection.pathogen, time_of_exposure, infected_agents)
                    #println("**agent $(source.id) infected agent $(target.id) with $pathogen_name")
                    infection_occurred =  true
                end
            end

        end
    end 

    return infection_occurred
end


# infect index case
function select_random_general_staff(agents::Agents, day::Int64)::Int64

    ids = collect(keys(agents.workers_G))
    
    index_case = 0 
    present = false
    n_max = 1000
    n = 0
    while !present
        index_case = sample(rng_infections, ids)
        present = (agents.workers_G[index_case].roster[day] == 1)
        n += 1
        if n > n_max
            println("trouble assigning random index worker index case, check rosters")
            return 0
        end
    end
    return index_case
end

function select_random_resident(agents)::Int64

    resident_ids = collect(keys(agents.residents))
    index_case_id = sample(rng_infections, resident_ids)
    return index_case_id 

end

function select_random_worker(agents::Agents, day::Int64)::Int64

    ids_g = collect(keys(agents.workers_G))
    ids_m = collect(keys(agents.workers_M))

    ids = vcat(ids_g, ids_m)
    
    index_case = 0 
    present = false
    n_max = 1000
    n = 0
    while !present
        index_case = sample(rng_infections, ids)
        present = (agents.All[index_case].roster[day] == 1)
        n += 1
        if n > n_max
            println("trouble assigning random index worker index case, check rosters")
            return 0
        end
    end
    return index_case
end

#NOTE: first we see if it's a worker, then we make sure they're present
# this ensures that p(index case is worker) = p(worker), and that 
# the definition of 'outbreak' doesn't change (i.e., the outbreak simulation)
# begins when the facility is exposed. 
function select_random_agent(agents::Agents, day::Int64)::Int64

    present = false
    ids = collect(keys(agents.All))
    index_case = sample(rng_infections, ids)
    if is_resident(agents.All[index_case])
        present = true 
    end
    n_max = 1000
    n = 0
    while !present
        index_case = sample(rng_infections, ids)
        present = (agents.All[index_case].roster[day] == 1)
        n += 1
        if n > n_max
            println("trouble assigning random index worker index case, check rosters")
            return 0
        end
    
    end
    return index_case
end


# update agent infections 
# NOTE: this function uses the global parameter dt. 
function update_infections!(agents::Agents, infected_agents::Dict{Int64, Float64})

    agents_fully_recovered = []
    for (id_i, t_inf) in infected_agents

        infections_to_remove = []

        # iterate infection foward in time 
        for (pathogen_name_j, infection_j) in agents.All[id_i].infections
            recovered = update_infection!(infection_j, dt)
            if recovered
                push!(infections_to_remove, pathogen_name_j)
            end
        end
        # remove any recovered infections 
        for pathogen_name in infections_to_remove
            #update immunity status for pathogen: 
            im = agents.All[id_i].immunity
            im.protection_Death[pathogen_name] = 1.0
            im.protection_Infection[pathogen_name] = 1.0
            im.protection_Symptoms[pathogen_name] = 1.0
            im.protection_Transmission[pathogen_name] = 1.0
            
            # HACK: for analysis of single outbreaks, implementing this as SIR: 
            # TODO: generalise the implementation for long-term multistrain dynamics. 


            delete!(agents.All[id_i].infections, pathogen_name)
            #println("agent $id_i recovered from $pathogen_name")

        end

        #check if agent is no longer infected: 
        if isempty(agents.All[id_i].infections)
            push!(agents_fully_recovered, id_i)
        end

    end
    # remove any infected agents who just recovered from all infections: 
    for key in agents_fully_recovered
        delete!(infected_agents, key)
    end
end



# sample infectious contacts (network-based transmission)
function compute_transmission!(all_transmissions::DataFrame, infected_agents:: Dict{Int64, Float64}, agents::Agents,
                               day_of_week::Int64, w_tot_d::Dict{Int64, Float64}, 
                               contact_rate::Float64, bkg_contact_rate::Float64, t::Float64)

    E_list_infectious_t = E_list()
    #iterate over infected agents and add edges to infectious E_list
    infected_resident_ids = Vector{Int64}() # for background contact sampling 
    for (id, t_inf) in infected_agents

        # if agent is removed, do not add their edges 
        if agents.All[id].t_removed < 0.0

            a = agents.All[id] 
            if haskey(a.contacts, day_of_week)
                if is_resident(a) # queuing for background contacts. 
                    push!(infected_resident_ids, a.id)
                end 
                add_source_edges_to_E_list!(E_list_infectious_t, a, a.contacts[day_of_week])
            end
        end
    end

    weights_infectious_edges_t = compile_weights( E_list_infectious_t )

    #sum weights of infected edges
    w_infected = sum_edge_weights_EList(E_list_infectious_t)
    prop_infected = w_infected / w_tot_d[day_of_week]

    # n infectious edges to sample: 

    # Poisson(net contact rate * prop_infected * dt) 
    infectious_contact_rate = contact_rate * prop_infected
    dist = Poisson(infectious_contact_rate)
    n_to_sample = rand(rng_infections, dist)

    edges_to_evaluate = sample_E_list(E_list_infectious_t, 
                                      n_to_sample, 
                                      weights_infectious_edges_t)
    
    add_background_contacts!(edges_to_evaluate, infected_resident_ids, agents, bkg_contact_rate )
    
    
                                      #= some debugging printouts
        #println("contact rate per step: $contact_rate_per_step")
        #println("infectious contact rate per step: $infectious_contact_rate")
        #println("p_infected at time $t : $prop_infected")
        #println("going to evaluate: $n_to_sample infectious edges")
    =#
    # compute pairwise transmission over edges: 
    # there is still the possibility of selecting 
    # edges between infected individuals, so we'll have to exclude those: 
    for e in edges_to_evaluate.edges

        source = agents.All[e.source_id]
        target = agents.All[e.target_id]


        if ( ACTIVE_OUTBREAK && OUTBREAK_CONTROL ) 
            transmission_occurred = transmit_infection_AOB!(source, target, t, infected_agents)
        else
            transmission_occurred = transmit_infection!(source, target, t, infected_agents)
        end
       
        if transmission_occurred
            if is_worker(target)
                if target.is_medical 
                    push!(all_transmissions, (source.id, target.id, t, 3))
                else 
                    push!(all_transmissions, (source.id, target.id, t, 2))
                end
            else
                push!(all_transmissions, (source.id, target.id, t, 1))
            end

        end
    end

end

function add_background_contacts!(edges_out::E_list, source_ids::Vector{Int64}, agents::Agents, bkg_contact_rate::Float64)

    dist = Poisson(bkg_contact_rate)
    resident_ids = collect(keys(agents.residents))

    for source_id in source_ids
        n_to_sample = rand(rng_infections, dist)
        if n_to_sample > 0
            n_sampled = 0
            while n_sampled < n_to_sample
                
                target_id = sample(rng_infections, resident_ids) # sampling background contacts from residents only
                if (target_id != source_id)
                    new_edge = edge_type(source_id, target_id, 1.0) #this will require memory allocation, may be faster to look up. 
                    push!(edges_out.edges, new_edge)
                    n_sampled += 1
                end
            end
        end
    end

end



# simulate testing and detection 

function test_agent(a::Agent_T)::Bool

    positive = false 

    for (pathogen_name, infection) in a.infections
        compute_test_sensitivity!(infection)
        if rand(rng_testing) < infection.test_sensitivity 
            positive = true
        end
    end

    return positive
end

function test_agents!(detections::DataFrame, p_test_per_day::Float64, 
                      agents::Agents, infected_agents::Dict{Int64, Float64}, day_of_week::Int64, t::Float64)

    # test for existing infections (RAT): 
    for (id, t_inf) in infected_agents
        a = agents.All[id] 

        # TODO implement detection upon symptom expression as well as testing 
        if haskey(a.contacts, day_of_week) # if they're present 
            if rand(rng_testing) < p_test_per_day
                tested_positive = test_agent(a)
                if tested_positive 
                    a.t_detected = t

                    type_index = 0

                    if is_worker(a)
                        if a.is_medical
                            type_index = 3
                        else
                            type_index = 2
                        end
                    else
                        type_index = 1
                    end

                    push!(detections, (a.id, t, type_index))
                end
            end
        end
    end
end

function test_workers!(detections::DataFrame, p_test_per_day::Float64, 
                       agents::Agents, infected_agents::Dict{Int64, Float64},
                       day_of_week::Int64, t::Float64, removed_agents::Dict{Int64, Float64})

    # test for existing infections (RAT): 
    for (id, t_inf) in infected_agents
        # check if they've already been detected: 
        

        a = agents.All[id] 
        if a.t_detected < 0.0
            if is_worker(a)
            # TODO implement detection upon symptom expression as well as testing 
                if (haskey(a.contacts, day_of_week) && (a.t_removed < 0.0 ) )# if they're present 
                    
                    p = p_test_per_day
                    if check_symptom_expression(a)
                        p = p_test_if_symptomatic
                    end
                    
                    if rand(rng_testing) < p
                        tested_positive = test_agent(a)
                        if tested_positive 
                            a.t_detected = t

                            #remove worker agent: 
                            remove_worker!(a, t, removed_agents)
                            

                            type_index = 0
                            if a.is_medical
                                type_index = 3
                            else
                                type_index = 2
                            end

                            push!(detections, (a.id, t, type_index))
                        end
                    end
                end
            end
        end
    end
end

function test_residents!(detections::DataFrame, p_test_per_day::Float64, 
                         agents::Agents, infected_agents::Dict{Int64, Float64}, 
                         day_of_week::Int64, t::Float64)

    # test for existing infections (RAT): 
    for (id, t_inf) in infected_agents
        a = agents.All[id] 
        if a.t_detected < 0.0 # initialised negative for this purpose 
            if is_resident(a)
            # TODO implement detection upon symptom expression as well as testing 
                if haskey(a.contacts, day_of_week) # if they're present 

                    p = p_test_per_day
                    if check_symptom_expression(a)
                        p = p_test_if_symptomatic
                    end
                    
                    if rand(rng_testing) < p

                        tested_positive = test_agent(a)
                        if tested_positive 
                            a.t_detected = t

                            type_index = 1

                            push!(detections, (a.id, t, type_index))
                        end
                    end
                end
            end
        end
    end
end

function check_symptom_expression(a::Agent_T)::Bool

    symptoms = false 
    for (name, im) in a.infections
        if im.expressing_symptoms
            symptoms = true
        end
    end

    return symptoms

end

# toggle active outbreak

function count_detections_last_7_days(detections::DataFrame, t)::Int64
   
    d = floor(t) 
    dm7 = d - 7

    detections_last_7 = detections[detections.time_detected .> dm7, :]

    return size(detections_last_7, 1)

end


# remove a worker : NOTE: don't need to modify w_tot_d IF we assume that removing a worker
# reduces the overall contact rate proportionally to w_i / w_tot (differences cancel)
# ALTERNATELY: if we re-distributed the contacts to other workers, the overall contact rate
# would not change and neither would the denominators in w_tot_d. 
# TODO: determine if partial redistribution of contacts would require altering denominators 
function remove_worker!(a::Agent_T, t_removed::Float64, removed_agents::Dict{Int64, Float64})
    # a removed worker can't infect or be infected. 
    if WORKER_CASE_ISOLATION
        a.t_removed = t_removed
        removed_agents[a.id] = t_removed
    end
end

function update_absentees!(agents::Agents, removal_period::Float64, removed_agents::Dict{Int64, Float64}, t::Float64 )

    n_absent = 0 

    agents_to_reinstate = []
    for (id_i, t_removed) in removed_agents
        if t - t_removed > removal_period
            a = agents.All[id_i]
            a.t_removed = -1.0 #resident
            push!(agents_to_reinstate, id_i)
        else
            n_absent += 1
        end
    end

    for id_i in agents_to_reinstate
        delete!(removed_agents, id_i)
    end

    #println("number of agents absent: $n_absent")

end
# reinstate a worker (modifies the contact weight denominators)