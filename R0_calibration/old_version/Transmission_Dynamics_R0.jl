

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
function select_random_general_staff(agents::Agents_T, day::Int64)::Int64

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

function select_random_worker(agents::Agents_T, day::Int64)::Int64

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
function select_random_agent(agents::Agents_T, day::Int64)::Int64

    ids = collect(keys(agents.All))
    index_case = sample(rng_infections, ids)

    # ensures we select a worker who is present,
    # while selecting workers with probability proportional 
    # to the worker fraction.
    if is_worker(agents.All[index_case])
        index_case = select_random_worker(agents, day)
    end

    return index_case
end


# update agent infections 
# NOTE: this function uses the global parameter dt. 
function update_infections!(agents::Agents_T, infected_agents::Dict{Int64, Float64})

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
function compute_transmission!(all_transmissions::DataFrame, infected_agents:: Dict{Int64, Float64}, agents::Agents_T,
                               day_of_week::Int64, w_tot_d::Dict{Int64, Float64}, 
                               contact_rate::Float64, bkg_contact_rate::Float64, bkg_contact_rate_iso::Float64, t::Float64)

    E_list_infectious_t = E_list()
    #iterate over infected agents and add edges to infectious E_list
    infected_resident_ids = Vector{Int64}() # for background contact sampling 
    for (id, t_inf) in infected_agents
        a = agents.All[id] 
        if haskey(a.contacts, day_of_week)
            if is_resident(a) # queuing for background contacts. 
                push!(infected_resident_ids, a.id) 
                #NOTE: 2022 09 19 nesting this under the isolation test means isolated residents were not added
                # this should now be fixed 
            end 

            active_contacts = Array{Contact_T, 1}()
            for c in a.contacts[day_of_week]
                push!(active_contacts, c) 
                # active contacts will include same-room resident contacts and any 
                # worker contacts where the target is not removed (i.e., furloughed)
                # NOTE: isolated residents still contribute to these contacts, 
                # their isolation status is taken into account w.r.t. background contacts
                # only. 
            end

            #if day_of_week == 5 && a.id == 155
            #    println("check here")
            #end

            add_source_edges_to_E_list!(E_list_infectious_t, a, active_contacts)#a.contacts[day_of_week])
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
    n_to_sample = rand(rng_contacts, dist)

    edges_to_evaluate = sample_E_list(E_list_infectious_t, 
                                      n_to_sample, 
                                      weights_infectious_edges_t)
    
    add_background_contacts!(edges_to_evaluate, infected_resident_ids, agents, bkg_contact_rate, bkg_contact_rate_iso )
    
    
                                      #= some debugging printouts
        #println("contact rate per step: $contact_rate_per_step")
        #println("infectious contact rate per step: $infectious_contact_rate")
        #println("p_infected at time $t : $prop_infected")
        #println("going to evaluate: $n_to_sample infectious edges")
    =#

    # 2022 09 19 testing: 
    all_edges_good = test_contacts(edges_to_evaluate, agents)
    if !all_edges_good
        println("WARNING: invalid infectious contacts found! Check network.")
    end


    # compute pairwise transmission over edges: 
    # there is still the possibility of selecting 
    # edges between infected individuals, so we'll have to exclude those: 
    for e in edges_to_evaluate.edges

        source = agents.All[e.source_id]
        target = agents.All[e.target_id]

        # NOTE: adding PPE flag to this condition (2022 09 23)
        if ( ACTIVE_OUTBREAK && OUTBREAK_CONTROL && PPE_AVAILABLE) 
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

function compute_transmission!(all_transmissions::DataFrame, infected_agents:: Dict{Int64, Float64}, agents::Agents_T,
                               day_of_week::Int64, w_tot_d::Dict{Int64, Float64}, 
                               contact_rate::Float64, bkg_contact_rate::Float64, bkg_contact_rate_iso::Float64, t::Float64)

    E_list_infectious_t = E_list()
    #iterate over infected agents and add edges to infectious E_list
    infected_resident_ids = Vector{Int64}() # for background contact sampling 
    for (id, t_inf) in infected_agents
        a = agents.All[id] 
        if haskey(a.contacts, day_of_week)
            if is_resident(a) # queuing for background contacts. 
                push!(infected_resident_ids, a.id) 
                #NOTE: 2022 09 19 nesting this under the isolation test means isolated residents were not added
                # this should now be fixed 
            end 

            active_contacts = Array{Contact_T, 1}()
            for c in a.contacts[day_of_week]
                push!(active_contacts, c) 
                # active contacts will include same-room resident contacts and any 
                # worker contacts where the target is not removed (i.e., furloughed)
                # NOTE: isolated residents still contribute to these contacts, 
                # their isolation status is taken into account w.r.t. background contacts
                # only. 
            end

            #if day_of_week == 5 && a.id == 155
            #    println("check here")
            #end

            add_source_edges_to_E_list!(E_list_infectious_t, a, active_contacts)#a.contacts[day_of_week])
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
    n_to_sample = rand(rng_contacts, dist)

    edges_to_evaluate = sample_E_list(E_list_infectious_t, 
                                      n_to_sample, 
                                      weights_infectious_edges_t)
    
    add_background_contacts!(edges_to_evaluate, infected_resident_ids, agents, bkg_contact_rate, bkg_contact_rate_iso )
    
    
                                      #= some debugging printouts
        #println("contact rate per step: $contact_rate_per_step")
        #println("infectious contact rate per step: $infectious_contact_rate")
        #println("p_infected at time $t : $prop_infected")
        #println("going to evaluate: $n_to_sample infectious edges")
    =#

    # 2022 09 19 testing: 
    all_edges_good = test_contacts(edges_to_evaluate, agents)
    if !all_edges_good
        println("WARNING: invalid infectious contacts found! Check network.")
    end


    # compute pairwise transmission over edges: 
    # there is still the possibility of selecting 
    # edges between infected individuals, so we'll have to exclude those: 
    for e in edges_to_evaluate.edges

        source = agents.All[e.source_id]
        target = agents.All[e.target_id]

        # NOTE: adding PPE flag to this condition (2022 09 23)
        if ( ACTIVE_OUTBREAK && OUTBREAK_CONTROL && PPE_AVAILABLE) 
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

function compute_transmission_R0!(index_id::Int64, all_transmissions::DataFrame, infected_agents:: Dict{Int64, Float64}, agents::Agents_T,
                               day_of_week::Int64, w_tot_d::Dict{Int64, Float64}, 
                               contact_rate::Float64, bkg_contact_rate::Float64, bkg_contact_rate_iso::Float64, t::Float64)

    E_list_infectious_t = E_list()
    #iterate over infected agents and add edges to infectious E_list
    infected_resident_ids = Vector{Int64}() # for background contact sampling 
        
    #for R0 estimates, we only examine the index case    
    id = index_id 
    if haskey(infected_agents, id) #if index case is still infected (termination flag is updated daily, so there are a few instances where index case is recovered. )

        a = agents.All[id] 
        if haskey(a.contacts, day_of_week)
            if is_resident(a) # queuing for background contacts. 
                push!(infected_resident_ids, a.id) 
                #NOTE: 2022 09 19 nesting this under the isolation test means isolated residents were not added
                # this should now be fixed 
            end 

            active_contacts = Array{Contact_T, 1}()
            for c in a.contacts[day_of_week]
                push!(active_contacts, c) 
                # active contacts will include same-room resident contacts and any 
                # worker contacts where the target is not removed (i.e., furloughed)
                # NOTE: isolated residents still contribute to these contacts, 
                # their isolation status is taken into account w.r.t. background contacts
                # only. 
            end

            #if day_of_week == 5 && a.id == 155
            #    println("check here")
            #end

            add_source_edges_to_E_list!(E_list_infectious_t, a, active_contacts)#a.contacts[day_of_week])
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
    n_to_sample = rand(rng_contacts, dist)

    edges_to_evaluate = sample_E_list(E_list_infectious_t, 
                                      n_to_sample, 
                                      weights_infectious_edges_t)
    
    add_background_contacts!(edges_to_evaluate, infected_resident_ids, agents, bkg_contact_rate, bkg_contact_rate_iso )
    
    
                                      #= some debugging printouts
        #println("contact rate per step: $contact_rate_per_step")
        #println("infectious contact rate per step: $infectious_contact_rate")
        #println("p_infected at time $t : $prop_infected")
        #println("going to evaluate: $n_to_sample infectious edges")
    =#

    # 2022 09 19 testing: 
    all_edges_good = test_contacts(edges_to_evaluate, agents)
    if !all_edges_good
        println("WARNING: invalid infectious contacts found! Check network.")
    end


    # compute pairwise transmission over edges: 
    # there is still the possibility of selecting 
    # edges between infected individuals, so we'll have to exclude those: 
    for e in edges_to_evaluate.edges

        source = agents.All[e.source_id]
        target = agents.All[e.target_id]

        # NOTE: adding PPE flag to this condition (2022 09 23)
        if ( ACTIVE_OUTBREAK && OUTBREAK_CONTROL && PPE_AVAILABLE) 
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


# testing utility: 
function test_contacts(edges_to_evaluate::E_list_T, agents::Agents_T)::Bool

    test_flag = true 

    for e in edges_to_evaluate.edges 

        s = agents.All[e.source_id] 
        t = agents.All[e.target_id]
        
        # tests: 
        #(1) is source infectious? 
        if !is_infected(s, "Omicron")
            println("agent $(s.id) is trying to infect agent $(t.id), but $(s.id) is not infected")
            test_flag = false 
        end

        #(2) is source an isolated worker? 
        if (is_worker(s) && is_isolated(s))
            println("agent $(s.id) is trying to infect agent $(t.id), but $(s.id) is an isolated worker")
            test_flag = false 
        end

        #(3) is target an isolated worker? 
        if (is_worker(t) && is_isolated(t))
            println("agent $(s.id) is trying to infect agent $(t.id), but $(t.id) is an isolated worker")
            test_flag = false 
        end

    end

    return test_flag 

end


function add_background_contacts!(edges_out::E_list_T, source_ids::Vector{Int64}, agents::Agents_T, bkg_contact_rate::Float64, bkg_contact_rate_iso::Float64)

    dist = Poisson(bkg_contact_rate)
    dist_iso = Poisson(bkg_contact_rate_iso)

    resident_ids = collect(keys(agents.residents))
    # sample should be weighted by isolation status 
    # NOTE: isolation reduces number of contacts made 
    # also reduces relative chance of being contacted
    iso_weights = Float64[]
    for r_id in resident_ids 
        # check isolation status 
        if !is_isolated(agents.residents[r_id])
            # resident is not isolated 
            if ACTIVE_OUTBREAK && RESIDENT_LOCKDOWN
            push!(iso_weights, 1.0 - resident_lockdown_efficacy)
            else 
                push!(iso_weights, 1.0)
            end
        else
            push!(iso_weights, 1.0 - resident_isolation_efficacy)
        end
    end

    for source_id in source_ids

        # not isolated 
        if !is_isolated(agents.All[source_id])
            n_to_sample = rand(rng_contacts, dist)
        else
            n_to_sample = rand(rng_contacts, dist_iso)
        end
        
        if n_to_sample > 0
            n_sampled = 0
            while n_sampled < n_to_sample
                
                target_id = sample(rng_contacts, resident_ids, Weights(iso_weights)) # sampling background contacts from residents only
                
                if (target_id != source_id)
                    new_edge = Edge(source_id, target_id, 1.0) #this will require memory allocation, may be faster to look up. 
                    push!(edges_out.edges, new_edge)
                    n_sampled += 1
                end
            end
        end
    end

end



