# Author: Cameron Zachreson
# Institution: The University of Melbourne
# Simulation code acompanying the manuscript entitled: 
# "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
# Date released: Dec. 18, 2023




global NETWORK_TEST = false 

# define the main function 

#### runs transmission dynamics with outbreak detection, until 
# a specified number of outbreaks have been recorded
# in the output linelist
function run_OB!(config::Setup_RACF.Config_T, 
                 pop::Setup_RACF.Pop_Input_T, 
                 facility::Facility_Structure.Facility_T,
                 n_tot::Int64,
                 n_outbreaks_flag::Bool,
                 n_instances_flag::Bool,
                 output_dirname::String) # ! because run can modify config. 
    
    #output linelist vectors 

    tot_infections = Vector{Int64}()
    tot_detections = Vector{Int64}()
    time_to_detection = Vector{Float64}()
    tot_infections_at_detection = Vector{Float64}()

    output_linelist = DataFrame(
        run_id = Int64[], 
        facility_id = Int64[],
        n_residents = Int64[],
        n_staff = Int64[],
        total_FTE = Float64[],
        index_case_type = String[],
        t_first_detection = Float64[],
        OB_declared = Bool[],
        t_OB_on = Float64[],
        t_OB_off = Float64[],
        I_cum_OB_on = Int64[],
        Det_cum_OB_on = Int64[],
        F_max = Int64[],
        Iso_max = Int64[],
        I_tot_staff = Int64[],
        Det_tot_staff = Int64[],
        I_tot_res = Int64[],
        Det_tot_res = Int64[],
        FTE_def_max = Float64[],
        sim_duration = Float64[]
    )


    #setup disease parameters: 
    diseases = Diseases_RACF.diseases()
    Diseases_RACF.set_disease_dict!(diseases, config)

    n = 0

    i = 0
    while (n < n_tot)
        #i in 1:n_runs

        i += 1
        config.active_outbreak = false

        if n_instances_flag
            n += 1
            println("\n****run:$i")
        end

        seed_offset = i # used for re-seeding rngs for each run, for comparability between each linelist element. 

        Setup_RACF.apply_seed_offset!(config, seed_offset) #should re-seed random number generators each with a new run-specific value. 
        # this will allow direct comparison of linelist entries from different scenarios on the same facility. 


        # fill any constant params into linelist: 
        push!(output_linelist.run_id, i)
        push!(output_linelist.facility_id, facility.id)
        

        #note the below are placeholders, that will be modified if the appropriate triggers are activated:
        push!(output_linelist.OB_declared, false)
        push!(output_linelist.t_OB_on, NaN)
        push!(output_linelist.t_OB_off, NaN)
        push!(output_linelist.I_cum_OB_on, -1) #note: no nan type for integers 
        push!(output_linelist.Det_cum_OB_on, -1)

            
        F_max = 0 # maximum number of concurrently furloughed staff 
        Iso_max = 0 # maximum number of isolated residents 
        FTE_def_max = 0 # maximum FTE deficit from furlough 

        ## ***** ##
        #initialise Agents 
        agents = Agents_RACF.Agents()


        ln_neut_dist_res = Normal(config.log_mu_res, config.log_sig_res)
        ln_neut_dist_wG = Normal(config.log_mu_wG, config.log_sig_wG)
        ln_neut_dist_WM = Normal(config.log_mu_wM, config.log_sig_wM)

        Agents_RACF.populate_workers_from_DataFrame_imDist!(agents.workers_G, pop.workers_G_str, ln_neut_dist_wG, config)#, N_lists.id_to_contacts)
        N_workers_G = length(agents.workers_G) #note size() not defined for dict
        #println("created $N_workers_G general staff members")

        #initialise medical staff
        Agents_RACF.populate_workers_from_DataFrame_imDist!(agents.workers_M, pop.workers_M_str, ln_neut_dist_WM, config)#, N_lists.id_to_contacts)
        N_workers_M = length(agents.workers_M) #note size() not defined for dict
        #println("created $N_workers_M medical staff members")

        #initialise residents
        Agents_RACF.populate_residents_from_DataFrame_imDist!(agents.residents, pop.residents_str, ln_neut_dist_res, config)#, N_lists.id_to_contacts)
        N_residents = length(agents.residents) #note size() not defined for dict
        #println("created $N_residents residents")


        agents.All = merge(agents.residents, 
                           agents.workers_G, 
                           agents.workers_M)

        # add agent numbers to output list 
        push!(output_linelist.n_residents, N_residents)
        push!(output_linelist.n_staff, (N_workers_G + N_workers_M))
        push!(output_linelist.total_FTE, Agents_RACF.compute_total_FTE(agents))

        #populate rooms (new on 2022 09 13, need it for re-distribution of labour during 
        # surge rostering from furloughs)
        rooms = Facility_Structure.Rooms()
        Facility_Structure.populate_Rooms_from_Agents!(rooms, agents)

        # assign neighbour lists from rooms:
        N_lists_new = Networks_RACF.N_list_temporal_multigraph() 
        n_days = 7 #hacked in.
        for d in 1:n_days
            N_lists_new.day_to_N_list[d] = Networks_RACF.N_list()
            N_list_d = N_lists_new.day_to_N_list[d]
            Facility_Structure.populate_N_lists_d_from_Rooms!(rooms, N_list_d, d)
        end

        # assign agent contacts from N lists constructed above from rooms 
        Agents_RACF.assign_contacts!(agents, N_lists_new) 
        #note - this will completely replace the existing contact list 

        # adjust pairwise weights based on model-specific bias 
        # NOTE: adjusting weights in agent neighbour lists 
        # also adjusts the weights in N_lists (referencing)
        Agents_RACF.pairwise_weights!(agents, config)

        if NETWORK_TEST
            # TEST: 2022 09 20
            # copy the network, and test against final network
            # after simulation is complete and all workers 
            # have returned from furlough
            E_list_initial =  Networks_RACF.E_list()
            for (id, a) in agents.All
                for (d, c) in a.contacts
                    Agents_RACF.add_source_edges_to_E_list!(E_list_initial, a, c )
                end
            end
        end

        #initialise test schedule 
        Outbreak_Response.initialise_baseline_testing_schedule!(agents, config)

        ##*****##

        #println("adjusted pairwise weights")

        edge_lists_all = Networks_RACF.E_list_temporal()
        # iterate through agents and contacts and generate 
        # full edge lists (all possible contacts)
        Agents_RACF.fill_E_lists_all!(edge_lists_all, agents)
        #println("created edge lists for each day")

        # compute weight totals for each day (denominators used to compute sampling rate)
        w_tot_d = Dict{Int64, Float64}() # day -> w_tot
        for (d, EL) in edge_lists_all.day_to_E_list
            #w_tot_d[d] = sum_edge_weights_EList(EL) # directed
            w_tot_d[d] = Networks_RACF.sum_edge_weights_EList(EL) / 2.0 #undirected
        end
        # NOTE: currently sum_edge_weights_EList() double-counts
        # all interactions because edges are accessed 
        # by searching neighbour lists, need the factor
        # of two in the contact rate calculation
        # because 'infectious' edges are only counted
        # from the infected source. 

        # make a dict for infected agents id -> time infected
        infected_agents = Dict{Int64, Float64}()
        removed_workers = Dict{Int64, Float64}()
        isolated_residents = Dict{Int64, Float64}()

        # output structures 
        all_transmissions = DataFrame(source = [], 
                                      target = [], 
                                      time = [], 
                                      type_of_agent = [])

        all_detections = DataFrame(id = [], 
                                   time_detected = [], 
                                   type_of_agent = [])
        

        # initialise transmission simulation

        # initialise time keepers 
        t_i = 0.0
        t_f = 120.0 # ~4 months, maximum simulation time. 
        termination_flag = false 
        step_final = convert(Int64, ceil((t_f - t_i)/config.dt ))

        # TODO: should set the system up to start on an arbitrary day. 
        day_of_week_i = convert(Int64, floor(mod(t_i, 7))) + 1 # first day. 

        #dt is defined in setup # dt = 0.1 (2022 08 17)

        #index case id
        #index_case_id = minimum(keys(agents.workers_G)) #first worker (id 101 in this example) 

        if config.resident_index_case && config.worker_index_case
            index_case_id = Transmission_Dynamics.select_random_agent(agents, 
                                                                      day_of_week_i,
                                                                      config)
        elseif config.resident_index_case
            index_case_id = Transmission_Dynamics.select_random_resident(agents,
                                                                         config)
        elseif config.worker_index_case
            #index_case_id = select_random_general_staff(agents, day_of_week_i)
            index_case_id = Transmission_Dynamics.select_random_worker(agents, 
                                                                       day_of_week_i,
                                                                       config)
        end

        if Agents_RACF.is_worker(agents.All[index_case_id])
            push!(output_linelist.index_case_type, "worker")
        else
            push!(output_linelist.index_case_type, "resident")
        end



        time_infected = t_i

        Transmission_Dynamics.infect_agent!(agents.All[index_case_id], 
                                            diseases.dict["Default"], 
                                            time_infected, 
                                            infected_agents,
                                            config)


        # define some parameters defining contact rates 
        #contact_rate_per_resident_per_day = 10.0 #NOTE: moved to setup_RACF.jl

        # 2022 09 28 contact rate per high needs resident vs rate per regular-needs resident 
        n_residents_high_needs = Agents_RACF.count_high_needs_residents(agents)
        n_residents_reg_needs = Agents_RACF.count_reg_needs_residents(agents)

        contact_rate_per_day = 
            convert(Float64, n_residents_high_needs) * 
            config.contact_rate_per_resident_per_day_high_needs
        # applying the same factor of 3 applied to the contact weights between workers and residents. 
        contact_rate_per_day += 
            convert(Float64, n_residents_reg_needs) * (config.contact_rate_per_resident_per_day)

        #contact_rate_per_day = convert(Float64, N_residents) * contact_rate_per_resident_per_day
        
        contact_rate_per_step = contact_rate_per_day * config.dt 

        # background contacts between residents (i.e., in shared meal spaces)
        # these are not included in the structured facility model
        # the implied structure is a fully-connected graph with uniform edge weight
        #bkg_contact_rate_per_resident_per_day = 5.0 #NOTE: moved to setup_RACF.jl
        bkg_contact_rate_per_resident_per_step = config.bkg_contact_rate_per_resident_per_day * config.dt



        # 2022 09 23 : including a delay between outbreak declaration and the implementation of measures 
        time_since_outbreak_declared = 0
        config.PPE_available = false #this gets toggled on after delay. 

        sim_duration_i = 0.0

        #for numeric times:
        small_num = 0.0000000001
        step = 0
        while !termination_flag
            step += 1
            #println(step)
            
            if step == step_final
                termination_flag = true
            end

            t = convert(Float64, step) * config.dt 
            day_of_week = convert(Int64, ceil(mod(t-small_num, 7))) # day is an integer 
            day = convert(Int64, ceil(t - small_num)) # first day. 


            #println("time is $t")
            #println("day is $day")

            t_last = convert(Float64, step-1) * config.dt 
            day_last = convert(Int64, ceil(t_last - small_num)) # day of previous step

            new_day = false
            if day_last < day
                new_day = true
                #println("new day, day of week: $day_of_week")
            end


            # iterate through infected agents and update
            # infection status of each.  
            Transmission_Dynamics.update_infections!(agents, infected_agents, config)

            #check for symptoms: 

            # test with RATs 
            if new_day 

                if isempty(infected_agents)
                    termination_flag = true
                end

                # iterate furlough periods for workers and re-instate 
                # NOTE: 2022 09 15, now adds any re-instated workers back to their 
                # assigned rooms. 
                room_ids_to_update = Set{Int64}()
                Outbreak_Response.update_absentees!(agents, 
                                                    config.removal_period, 
                                                    removed_workers, 
                                                    t, 
                                                    rooms, 
                                                    room_ids_to_update)
                # room_ids_to_update will now include all rooms where workers have been
                # reinstated from furlough 

                #iterate isolation periods for residents and re-instate
                Outbreak_Response.update_isolated_residents!(agents, 
                                                             config.removal_period, 
                                                             isolated_residents, 
                                                             t)
                
                F_t = length(removed_workers)
                if F_t > F_max
                    F_max = F_t # maximum number of concurrently furloughed staff 
                end

                Iso_t = length(isolated_residents)
                if Iso_t > Iso_max
                    Iso_max = Iso_t # maximum number of isolated residents 
                end
                
                FTE_def_t = Agents_RACF.compute_FTE_deficit(agents, removed_workers)
                if FTE_def_t > FTE_def_max
                    FTE_def_max = FTE_def_t # maximum FTE deficit from furlough 
                end
                

                worker_ids_to_remove = Array{Int64, 1}()
                resident_ids_to_isolate = Array{Int64, 1}() 

                # check to see whether to implement infection control after delay: 
                if time_since_outbreak_declared >= config.delay_infection_control
                    config.PPE_available = true #toggled on after delay. #see trigger for transmit_infection_AOB in #transmission_dynamics.jl  
                end

                # update delay clock for introduction of infection control measures: 
                if config.active_outbreak
                    time_since_outbreak_declared += 1 # integer increase b/c it's in the new_day scope 
                end

                


                # NOTE: Do we assume RATs are available immediately after outbreak declaration, 
                # or should delay be applied to this as well? 
                # For now, I'll assume tests are immediately available, but PPE stockpile is not. 
                if (config.active_outbreak && config.outbreak_control)
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    Outbreak_Response.test_workers!(all_detections, 
                                                    config.p_test_per_day_workers_outbreak, 
                                                    agents, 
                                                    infected_agents, 
                                                    day_of_week, 
                                                    t, 
                                                    worker_ids_to_remove, 
                                                    config)
                    # test residents 
                    Outbreak_Response.test_residents!(all_detections, 
                                                      config.p_test_per_day_residents_outbreak, 
                                                      agents, 
                                                      infected_agents, 
                                                      day_of_week, 
                                                      t, 
                                                      resident_ids_to_isolate, 
                                                      config)
                else
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    Outbreak_Response.test_workers!(all_detections, 
                                                    config.p_test_per_day_workers_baseline, 
                                                    agents, 
                                                    infected_agents, 
                                                    day_of_week, 
                                                    t, 
                                                    worker_ids_to_remove, 
                                                    config)
                    # test residents 
                    Outbreak_Response.test_residents!(all_detections, 
                                                      config.p_test_per_day_residents_baseline, 
                                                      agents, 
                                                      infected_agents, 
                                                      day_of_week, 
                                                      t, 
                                                      resident_ids_to_isolate, 
                                                      config)
                end

                # remove workers from rooms and from N_lists if they tested positive 
                if config.worker_case_isolation
                    Outbreak_Response.remove_workers!(agents, 
                                                      worker_ids_to_remove, 
                                                      N_lists_new, 
                                                      removed_workers, 
                                                      t, 
                                                      rooms, 
                                                      room_ids_to_update)
                end

                # isolate residents who tested positive: 
                if config.resident_case_isolation
                    Outbreak_Response.isolate_residents!(agents, 
                                                         resident_ids_to_isolate, 
                                                         isolated_residents, 
                                                         t)
                end


                # update network based on current room assignments 
                # NOTE: this needs to be optimised to update only those agents who are affected
                # by fuloughs 
                
                # NOTE: all days must be updated because the set room_ids_to_update refreshes each new day 
                # but accounts for effects on the entire roster. 

            
                if !isempty(room_ids_to_update)
                    #println("attempting to re-wire")
                    for (d, N_list_d) in N_lists_new.day_to_N_list
                        # note, out-edges from removed workers have been removed in remove_workers!()
                        # the below removes the in-edges to removed workers 
                        # (and will perform any other updates based on room assignment changes in surge roster (for future version))
                        # this also updates out-edges for any workers who are returning from furlough. 
                        Facility_Structure.update_N_lists_d_from_Rooms!(rooms, 
                                                                        room_ids_to_update, 
                                                                        N_list_d, 
                                                                        d)
                    end
                end
                # assign agent contacts from N lists constructed above from rooms 
                # only doing this for the current day, as the future days will probably change again 
                # before the contacts are needed. 
                Agents_RACF.assign_todays_contacts!(agents, 
                                                    N_lists_new, 
                                                    day_of_week) 
                #note - this will completely replace the existing contact list for today
        
                # adjust pairwise weights based on model-specific bias 
                Agents_RACF.pairwise_weights_d!(agents, 
                                                day_of_week,
                                                config) 
                # for the daily network update, this should be 
                # incorporated into assign_todays_contacts!()

                
                # count detections from last 7 days, if > 0, we have an active outbreak, otherwise we don't_detected
                #println("\n****\nt=$t: new cases last 7 days: $n_cases_last_7")

                # declare active outbreak, or delcare active outbreak over. 
                if !(config.active_outbreak)
                    OB_flag = Outbreak_Response.declare_outbreak(all_detections, t)
                    if OB_flag
                        config.active_outbreak = true
                        if n_outbreaks_flag
                            n += 1
                            println("\n****run:$i")
                            println("active outbreak declared at time $t")
                            println("there are: $(size(all_transmissions,1)) infections ")
                            println("total outbreaks: $n")
                        end
                        output_linelist.OB_declared[i] = true 
                        output_linelist.t_OB_on[i] = t
                        output_linelist.I_cum_OB_on[i] = size(all_transmissions,1) + 1 # add one to include index case 
                        output_linelist.Det_cum_OB_on[i] = size(all_detections, 1) 
                    end
                else
                    OB_over_flag = Outbreak_Response.outbreak_over(all_detections, t)
                    if OB_over_flag
                        config.active_outbreak = false
                        println("active outbreak over at time $t")
                        output_linelist.t_OB_off[i] = t
                        termination_flag = true
                    end
                end

            end 
            
            # this now implements PPE when outbreak is active. 
            #modified intputs: all_transmissions, infected_agents, agents

            bkg_contact_rate = bkg_contact_rate_per_resident_per_step
            bkg_contact_rate_iso = bkg_contact_rate_per_resident_per_step * (1.0 - config.resident_isolation_efficacy)

            if config.active_outbreak && config.outbreak_control && config.resident_lockdown
                bkg_contact_rate *= (1.0 - config.resident_lockdown_efficacy)
            end

            Transmission_Dynamics.compute_transmission!(all_transmissions, 
                                                        infected_agents, 
                                                        agents,
                                                        day_of_week, 
                                                        w_tot_d, 
                                                        contact_rate_per_step, 
                                                        bkg_contact_rate, 
                                                        bkg_contact_rate_iso, 
                                                        t, 
                                                        config)

            if termination_flag
                sim_duration_i = t
            end
        end

        # 2022 09 20
        #TEST: build final network and compare to original one
        # (should be the same if all workers have been reinstated by tf)
        if NETWORK_TEST
            E_list_final =  Networks_RACF.E_list()
            for (id, a) in agents.All
                for (d, c) in a.contacts
                    Agents_RACF.add_source_edges_to_E_list!(E_list_final, a, c )
                end
            end

            test_val_e1_e2 = Networks_RACF.compare_E_lists(E_list_initial, E_list_final)
            if test_val_e1_e2 
                println("final network same as initial network e_i -> e_f")
            else
                println("final network different from initial network e_i -> e_f")
            end
            test_val_e2_e1 = Networks_RACF.compare_E_lists(E_list_final, E_list_initial)
            if test_val_e2_e1 
                println("final network same as initial network e_f -> e_i")
            else
                println("final network different from initial network e_f -> e_i")
            end
        end

        if config.write_transmission_tree_flag && output_linelist.OB_declared[i]
            rand_label = now()
            rand_label = replace(string(rand_label), ":"=> "_")
            rand_label = replace(string(rand_label), "."=> "_")
            rand_label = replace(string(rand_label), "-"=> "_")
            if config.outbreak_control
                ob_lab = "ICAOB_"
            else
                ob_lab = "unmitigated_"
            end
            
            output_dir_run = "$output_dirname\\run_$(i)"
            if !ispath(output_dir_run)
                mkpath(output_dir_run)
            end

            CSV.write("$(output_dir_run)\\$(ob_lab)detected_cases_$rand_label.csv", all_detections)
            CSV.write("$(output_dir_run)\\$(ob_lab)transmission_tree_$rand_label.csv", all_transmissions)
        end

        push!(output_linelist.F_max, F_max)# maximum number of concurrently furloughed staff 
        push!(output_linelist.Iso_max, Iso_max) # maximum number of isolated residents 
        push!(output_linelist.FTE_def_max, FTE_def_max) # maximum FTE deficit from furlough 


        push!(output_linelist.sim_duration, sim_duration_i)
        # count total infections and detections in staff and residents 

        I_tot_staff = Agents_RACF.count_worker_infections(agents, all_transmissions)
        I_tot_res = Agents_RACF.count_resident_infections(agents, all_transmissions)
        Det_tot_staff = Agents_RACF.count_worker_detections(agents, all_detections)
        Det_tot_res = Agents_RACF.count_resident_detections(agents, all_detections)

        if Agents_RACF.is_worker(agents.All[index_case_id])
            I_tot_staff += 1
        else 
            I_tot_res += 1
        end

        push!(output_linelist.I_tot_staff, I_tot_staff) 
        push!(output_linelist.Det_tot_staff, Det_tot_staff) 
        push!(output_linelist.I_tot_res, I_tot_res) 
        push!(output_linelist.Det_tot_res, Det_tot_res) 



        if isempty(all_transmissions)
            total_infections = 1
        else
            total_infections = size(all_transmissions, 1) + 1 # + 1 is for index case. 
        end

        #println("total infections for run: $total_infections")

        push!(tot_infections, total_infections)
        
        if !isempty(all_detections)   
            t_first_detection = minimum(all_detections.time_detected)
            total_detections = size(all_detections, 1)
            

            total_infections_at_detection = size(all_transmissions[all_transmissions.time .< t_first_detection, :], 1) + 1

        else
            total_detections = 0
            t_first_detection = NaN
            total_infections_at_detection = NaN
            
        end
        push!(output_linelist.t_first_detection, t_first_detection)
        push!(time_to_detection, t_first_detection)
        push!(tot_detections, total_detections)
        push!(tot_infections_at_detection, total_infections_at_detection)
    
    end

    println("avg. total infections: $(mean(tot_infections))")
    println("avg. total detected cases: $(mean(tot_detections))")
    println("avg. t first detection: $(mean(time_to_detection[.!isnan.(time_to_detection)]))")
    println("avg. infections at first detection: $(mean(tot_infections_at_detection[.!isnan.(tot_infections_at_detection)]))")

    output_fname = "$(output_dirname)\\output_summary.csv"
    output_linelist_fname = "$(output_dirname)\\output_linelist.csv"

    output = DataFrame(I_tot = tot_infections, 
                    I_detected = tot_detections, 
                    t_first_detection = time_to_detection, 
                    I_t_detection = tot_infections_at_detection)

    CSV.write(output_fname, output)
    CSV.write(output_linelist_fname, output_linelist)



end

function main_OB()

    # loop through exemplars and run the simulations: 

    output_dir_L1 = "" #placeholder

    # RANGE OF CONTROL PARAMETERS: 

    test_configurations = String["asymp_testing",
                                 "asymp_testing_OB_only",
                                 "no_asymp_testing", 
                                 "unmitigated"]
    
    lockdown_compliance = Float64[0.9, 0.75, 0.5, 0.25, 0.0]
    
    n_outbreaks_tot = 1000 # how many 'declared' outbreaks to simulate before terminating each run loop
    # for nice distributions, 1000 is a good number (takes about 1hr per sceneario)
    n_instances_tot = 10000# how man 'instances' to run if detection is turned off. 
    n = 0

    count_outbreaks_flag = true #if we're counting declared outbreaks
    count_instances_flag = false #if we're counting instances 
    
    # facility data for initialisation of ABM population: 
    data_dir_L1 = pwd() 
    
    data_dir_L2 = "input\\hypotheticals"
    
    # list of facilities
    facility_list_dir = "$data_dir_L1\\$data_dir_L2"
    facility_list_fname = "$facility_list_dir\\hypothetical_facility_characteristics.csv"
    facility_list = DataFrame(CSV.File(facility_list_fname, delim = ","))

    #facility 3 - no shared rooms, no high-needs residents
    facility_index = 3 
    facility_label = "facID_$(facility_list.service_id[facility_index])_hyp"
    data_dirname = "$(data_dir_L1)\\$(data_dir_L2)\\$facility_label"#"$(data_dir_L1)$(data_dir_L2)\\$group_label\\$OB_label"
    
    if !ispath(data_dirname)
        println("terminating simulation: facility info not found - check input data")
        return
    end

    output_dir_L1 = pwd() * "\\output_v9_OB_test\\$facility_label"
    if !ispath(output_dir_L1)
        mkpath(output_dir_L1)
    end

    for (TC_i) in test_configurations

        test_config_label = TC_i

        if TC_i == "unmitigated"
            lockdown_compliance_i = [0.0]
        else
            lockdown_compliance_i = lockdown_compliance
        end

    for LD_j in lockdown_compliance_i  

        LD_str = replace("$LD_j", "."=> "p")
        LD_label = "lockdown_compliance_$LD_str"

        

            

            output_dir_L2 = "\\$test_config_label\\$LD_label"
            
            output_dir_fac = "$(output_dir_L1)\\$(output_dir_L2)"
            if !ispath(output_dir_fac)
                mkpath(output_dir_fac)
            end

            config_run = Setup_RACF.run_configuration()
            Setup_RACF.setup_run_default!(config_run, data_dirname)

            # set config parameters
            Setup_RACF.set_immunity_dist!(config_run)
            
            #testing parameters: 

            if TC_i == "no_asymp_testing"
                config_run.p_test_per_day_workers_baseline = 0.0
                config_run.p_test_per_day_residents_baseline = 0.0
                config_run.p_test_per_day_workers_outbreak = 0.0
                config_run.p_test_per_day_residents_outbreak = 0.0
                config_run.p_test_if_symptomatic = 1.0

                n = n_outbreaks_tot
                count_outbreaks_flag = true #if we're counting declared outbreaks
                count_instances_flag = false #if we're counting instances 

            end

            if TC_i == "asymp_testing_OB_only"
                config_run.p_test_per_day_workers_baseline = 0.0
                config_run.p_test_per_day_residents_baseline = 0.0
                config_run.p_test_per_day_workers_outbreak = 1.0
                config_run.p_test_per_day_residents_outbreak = 1.0
                config_run.p_test_if_symptomatic = 1.0

                n = n_outbreaks_tot
                count_outbreaks_flag = true #if we're counting declared outbreaks
                count_instances_flag = false #if we're counting instances 
            end

            if TC_i == "asymp_testing"
                config_run.p_test_per_day_workers_baseline = 1.0
                config_run.p_test_per_day_residents_baseline = 1.0
                config_run.p_test_per_day_workers_outbreak = 1.0
                config_run.p_test_per_day_residents_outbreak = 1.0
                config_run.p_test_if_symptomatic = 1.0

                n = n_outbreaks_tot
                count_outbreaks_flag = true #if we're counting declared outbreaks
                count_instances_flag = false #if we're counting instances 
            end

            if TC_i == "unmitigated"
                config_run.p_test_per_day_workers_baseline = 0.0
                config_run.p_test_per_day_residents_baseline = 0.0
                config_run.p_test_per_day_workers_outbreak = 0.0
                config_run.p_test_per_day_residents_outbreak = 0.0
                config_run.p_test_if_symptomatic = 0.0

                n = n_instances_tot
                count_outbreaks_flag = false #if we're counting declared outbreaks
                count_instances_flag = true #if we're counting instances 
            end

            #distancing/lockdown parameters: 

            config_run.resident_lockdown_efficacy = LD_j

            config_run.resident_lockdown = true 
            config_run.worker_case_isolation = true 
            config_run.resident_case_isolation = true
            config_run.resident_isolation_efficacy = 0.9 # imperfect compliance with resident case isolation 
            # NOTE: case isolation efficacy must be leq (<=) lockdown efficacy or background contact
            # assignment may fail TODO: assert this constraint or alter implementation to accept deviations
            config_run.removal_period = 7 #7-day furlough for positive workers

            #update the config parameters (ensuring any interdependent parameters are changed): 

            Setup_RACF.update_config!(config_run)

            #write record of run configuration 
            Setup_RACF.write_config_details(config_run, output_dir_fac)

            pop_run = Setup_RACF.population_input()
            Setup_RACF.read_in_population_data!(pop_run, config_run)

            facility = Facility_Structure.facility()
            facility.id = facility_list.service_id[facility_index]

            run_OB!(config_run, pop_run, facility, n, count_outbreaks_flag, count_instances_flag, output_dir_fac)
        end

    end
end




