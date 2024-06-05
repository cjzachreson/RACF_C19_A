# Author: Cameron Zachreson
# Institution: The University of Melbourne
# Simulation code acompanying the manuscript entitled: 
# "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
# Date released: Dec. 18, 2023



global NETWORK_TEST = false 

# define the main function 
###### facility-specific R0 calibration: 
function run_R0!(config::Setup_RACF.Config_T, 
              pop::Setup_RACF.Pop_Input_T, 
              facility::Facility_Structure.Facility_T,
              n_instances_tot::Int64,
              output_dirname::String)::Float64 # ! because run can modify config. 
    
    #output linelist vectors 

    tot_transmissions = Vector{Int64}()


    output_linelist = DataFrame(
        run_id = Int64[], 
        facility_id = Int64[],
        n_residents = Int64[],
        n_staff = Int64[],
        total_FTE = Float64[],
        index_case_type = String[],
        I_tot_staff = Int64[],
        I_tot_res = Int64[],
        index_case_id = Int64[],
        index_case_beta_max = Float64[],
        t_index_case_recovery = Float64[]
    )


    #setup disease parameters: 
    diseases = Diseases_RACF.diseases()
    Diseases_RACF.set_disease_dict!(diseases, config)


    # initialise agents here, instead of within the instance loop. 
    # note, agent states are reset in the instance loop. 

       ## ***** ##
    #initialise Agents 
    agents = Agents_RACF.Agents()

    ln_neut_dist_res = Normal(config.log_mu_res, config.log_sig_res)
    ln_neut_dist_wG = Normal(config.log_mu_wG, config.log_sig_wG)
    ln_neut_dist_wM = Normal(config.log_mu_wM, config.log_sig_wM)

    Agents_RACF.populate_workers_from_DataFrame_imDist!(agents.workers_G, pop.workers_G_str, ln_neut_dist_wG, config)#, N_lists.id_to_contacts)
    N_workers_G = length(agents.workers_G) #note size() not defined for dict
    #println("created $N_workers_G general staff members")

    #initialise medical staff
    Agents_RACF.populate_workers_from_DataFrame_imDist!(agents.workers_M, pop.workers_M_str, ln_neut_dist_wM, config)#, N_lists.id_to_contacts)
    N_workers_M = length(agents.workers_M) #note size() not defined for dict
    #println("created $N_workers_M medical staff members")

    #initialise residents
    Agents_RACF.populate_residents_from_DataFrame_imDist!(agents.residents, pop.residents_str, ln_neut_dist_res, config)#, N_lists.id_to_contacts)
    N_residents = length(agents.residents) #note size() not defined for dict
    #println("created $N_residents residents")

    agents.All = merge(agents.residents, 
                        agents.workers_G, 
                        agents.workers_M)

    #***** set network *****
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

    n_instances = 0
    i = 0
    while n_instances < n_instances_tot
        #i in 1:n_runs

        n_instances += 1

        # reset agent states (but do not re-initialise)
        Agents_RACF.reset_agent_states!(agents.residents, ln_neut_dist_res, config)
        Agents_RACF.reset_agent_states!(agents.workers_M, ln_neut_dist_wM, config)
        Agents_RACF.reset_agent_states!(agents.workers_G, ln_neut_dist_wG, config)
        #println("initialisd agent infection status for run $n_instances")

        config.active_outbreak = false

        i += 1
        #println("\n****run:$i")

        seed_offset = i # used for re-seeding rngs for each run, for comparability between each linelist element. 

        Setup_RACF.apply_seed_offset!(config, seed_offset) #should re-seed random number generators each with a new run-specific value. 
        # this will allow direct comparison of linelist entries from different scenarios on the same facility. 


        # fill any constant params into linelist: 
        push!(output_linelist.run_id, i)
        push!(output_linelist.facility_id, facility.id)
        

        # add agent numbers to output list 
        push!(output_linelist.n_residents, N_residents)
        push!(output_linelist.n_staff, (N_workers_G + N_workers_M))
        push!(output_linelist.total_FTE, Agents_RACF.compute_total_FTE(agents))


        #initialise test schedule 
        # note: test days are randomised within the instance loop 
        # though this is not relevant for R0... 
        # considering this implementation for the outbreak 
        # simulator. 
        Outbreak_Response.initialise_baseline_testing_schedule!(agents, config)

        ##*****##

   

        # make a dict for infected agents id -> time infected
        infected_agents = Dict{Int64, Float64}()

        # output structures 
        all_transmissions = DataFrame(source = [], 
                                      target = [], 
                                      time = [], 
                                      type_of_agent = [])

        
        # initialise transmission simulation

        # initialise time keepers 
        t_i = 0.0
        t_f = 90.0 # three months, maximum simulation time. 
        termination_flag = false 
        step_final = convert(Int64, ceil((t_f - t_i)/config.dt ))


        #dt is defined in setup # dt = 0.1 (2022 08 17)
        #index case id
        #note this returns a vector
        index_case_id = Agents_RACF.select_random_agents_by_weight(agents, config, 1)
        #index_case_id = Agents_RACF.select_random_agents_uniform(agents, config, 1)
        index_case_id = index_case_id[1]

        #println("\n****run:$i, index case: $index_case_id")

        push!(output_linelist.index_case_id, index_case_id)
        

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

        push!(output_linelist.index_case_beta_max, 
              agents.All[index_case_id].infections["Default"].beta_max)


        # define some parameters defining contact rates 
        #contact_rate_per_resident_per_day = 10.0 #NOTE: moved to setup_RACF.jl

        # 2022 09 28 contact rate per high needs resident vs rate per regular-needs resident 
        n_residents_high_needs = Agents_RACF.count_high_needs_residents(agents)
        n_residents_reg_needs = Agents_RACF.count_reg_needs_residents(agents)

        contact_rate_per_day = 
            convert(Float64, n_residents_high_needs) * 
            config.contact_rate_per_resident_per_day_high_needs

        contact_rate_per_day += 
            convert(Float64, n_residents_reg_needs) * (config.contact_rate_per_resident_per_day)
        
        contact_rate_per_step = contact_rate_per_day * config.dt 

        # background contacts between residents (i.e., in shared meal spaces)
        # these are not included in the structured facility model
        # the implied structure is a fully-connected graph with uniform edge weight
        #bkg_contact_rate_per_resident_per_day = 5.0 #NOTE: moved to setup_RACF.jl
        bkg_contact_rate_per_resident_per_step = config.bkg_contact_rate_per_resident_per_day * config.dt


        #for numeric times:
        small_num = 0.0000000001
        step = 0
        while !termination_flag
            step += 1
            
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

            if !Transmission_Dynamics.is_infected(agents.All[index_case_id], "Default")
                termination_flag = true 
                push!(output_linelist.t_index_case_recovery, t)
            end


            if new_day 

                if isempty(infected_agents)
                    termination_flag = true
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

            end 
            

            bkg_contact_rate = bkg_contact_rate_per_resident_per_step


            Transmission_Dynamics.compute_transmission_R0!(all_transmissions, 
                                                           infected_agents, 
                                                           agents,
                                                           day_of_week, 
                                                           w_tot_d, 
                                                           contact_rate_per_step, 
                                                           bkg_contact_rate, 
                                                           t, 
                                                           config,
                                                           index_case_id)

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

        if config.write_transmission_tree_flag
            rand_label = now()
            rand_label = replace(string(rand_label), ":"=> "_")
            rand_label = replace(string(rand_label), "."=> "_")
            rand_label = replace(string(rand_label), "-"=> "_")

            ob_lab = "unmitigated_"

            
            output_dir_run = "$output_dirname\\run_$(i)"
            if !ispath(output_dir_run)
                mkpath(output_dir_run)
            end
            CSV.write("$(output_dir_run)\\$(ob_lab)transmission_tree_$rand_label.csv", all_transmissions)
        end

        # count total infections and detections in staff and residents 

        I_tot_staff = Agents_RACF.count_worker_infections(agents, all_transmissions)
        I_tot_res = Agents_RACF.count_resident_infections(agents, all_transmissions)

        if Agents_RACF.is_worker(agents.All[index_case_id])
            I_tot_staff += 1
        else 
            I_tot_res += 1
        end

        push!(output_linelist.I_tot_staff, I_tot_staff) 
        push!(output_linelist.I_tot_res, I_tot_res) 

        push!(tot_transmissions, size(all_transmissions, 1))
        
    end

    println("avg. secondary cases: $(mean(tot_transmissions))")

    output_fname = joinpath("$(output_dirname)", "output_summary.csv")
    output_linelist_fname = joinpath("$(output_dirname)", "output_linelist.csv")

    output = DataFrame(secondary_cases_tot = tot_transmissions)

    CSV.write(output_fname, output)
    CSV.write(output_linelist_fname, output_linelist)

    return mean(tot_transmissions)


end

function main_R0()

    # loop through exemplars and run the simulations: 


    output_dir_L1 = "" #placeholder

    # DEFINE RANGE OF CONTROL PARAMETERS: 
    
    #vaccine-acquired immunity: 
    immunity_states = Bool[false]
    transmission_scalers = collect(0.005:0.005:0.5)
    
    n_instances_tot = 1000 # how many instances of the primary case simulation
    n_label = "n_$n_instances_tot"
    # for nice distributions, 1000 is a good number (takes about 1hr per sceneario)
    
    # 2022 09 19 : tests 
    # (2) check that no agents who are removed are tested [x]
    # iterate through outbreak-specific input files: 
    
    data_dir_L1 = pwd() #C:\\Users\\czachreson\\Desktop\\policy_work\\RACF_2022\\RACF_code\\population_generator\\2022_09_02"
    
    data_dir_L2 = joinpath("input", "hypotheticals")#"\\network_generator_step_4\\agents_OB_model_ready_exemplar\\2022_09_02"
    
    # list of facilitis to iterate through
    fac_list_dir = joinpath("$data_dir_L1", "$data_dir_L2")
    fac_list_fname = joinpath("$fac_list_dir", "hypothetical_facility_characteristics_revision.csv")
    fac_list = DataFrame(CSV.File(fac_list_fname, delim = ","))
    
    n_populations = size(fac_list, 1)


    for im in immunity_states 

        if im
            immunity_label = "immunity_on"
        else
            immunity_label = "immunity_off"
        end

        for i in 1:1 #facility 1 has half as many staff as facility 3 in the original list 

            fac_i = i 
            fac_label = "facID_$(fac_list.service_id[fac_i])_hyp"

            output_dir_L1 = joinpath(pwd(), 
                                     "output_v9_R0_test_rev_1", 
                                     "$immunity_label",
                                     "$fac_label", 
                                     "$n_label")

            if !ispath(output_dir_L1)
                mkpath(output_dir_L1)
            end

            secondary_cases = zeros(size(transmission_scalers))
            it = 0
            for tscale_j in transmission_scalers
                it += 1

                data_dirname = joinpath("$(data_dir_L1)", "$(data_dir_L2)", "$fac_label")

                if !ispath(data_dirname)
                    continue 
                end

                tscale_str = replace("$tscale_j", "."=> "p")
                tscale_label = "tscale_$tscale_str"



                output_dir_L2 = "$tscale_label"
                output_dir_fac = joinpath("$(output_dir_L1)", "$(output_dir_L2)")
                if !ispath(output_dir_fac)
                    mkpath(output_dir_fac)
                end

                config_run = Setup_RACF.run_configuration()
                Setup_RACF.setup_run_default!(config_run, data_dirname)

                # modify any of the default config
                # parameters that are adjusted by the run loop
                # globals.

                config_run.uniform_immunity = false
                
                Setup_RACF.set_immunity_dist!(config_run)
                
                config_run.transmission_scaler = tscale_j

                config_run.immunity = im

                # turn off case detection 
                config_run.p_test_per_day_residents_baseline = 0.0
                config_run.p_test_per_day_workers_baseline = 0.0
                config_run.p_test_if_symptomatic = 0.0

                #set contact rates: 
                config_run.contact_rate_per_resident_per_day = 3.0
                config_run.bkg_contact_rate_per_resident_per_day = 3.0
                config_run.contact_rate_per_resident_per_day_high_needs = 0.0 #no high needs residents anyway

                #update the config parameters (ensuring any interdependent parameters are changed): 

                Setup_RACF.update_config!(config_run)

                #write record of run configuration 
                Setup_RACF.write_config_details(config_run, output_dir_fac)

                pop_run = Setup_RACF.population_input()
                Setup_RACF.read_in_population_data!(pop_run, config_run)

                facility = Facility_Structure.facility()
                facility.id = fac_list.service_id[fac_i]

                println("*****running R0 simulator for <beta_max> = $tscale_j*****")
                
                R0_est = run_R0!(config_run, pop_run, facility, n_instances_tot, output_dir_fac)
                println("R0_est = $R0_est")
                println("**********")

                secondary_cases[it] = R0_est

            end

            output_main_R0 = DataFrame(beta = transmission_scalers, 
                                       R_est = secondary_cases)

            output_fname = joinpath("$output_dir_L1", "beta_vs_R_est.csv")

            CSV.write(output_fname, output_main_R0)

        end

    end

end



