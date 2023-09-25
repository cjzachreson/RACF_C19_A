# 2022 09 26 checking secondary case distribution 
# for index cases. 

#NOTE: initially, this test will assume no outbreak response 
# i.e., no case isolation etc. 

# this is switched off in setup_R0_RACF.jl 

# transmission from non-index cases is turned off. 
# this is done by ensuring the 'infected agents' 
# dict passed to compute_transmission!()
# only contains the index case. 



include("./header_R0_RACF.jl")

# DEFINE RANGE OF CONTROL PARAMETERS: 
# delay values in days 
global delay_vals = Int64[0]#, 6]#Int64[0, 2, 4, 6]
# compliance with scheduled testing (workers only)
global test_compliance = Float64[1.0]#, 0.75, 0.5, 0.25, 0.0]

global n_runs = 1000 # how many importations to simulate for the R0 estimate. 
## 
global delay_val_i = 0 # placeholder these get re-valued on loop 
global test_compliance_j = 1.0 # these get re-valued on loop 
# NOTE: test compliance parameter is implemented in the setup_RACF.jl script
global PPE_AVAILABLE = false #this gets toggled on after delay. 

#NOTE: moving these flags from setup to preamble. 
global IMMUNITY = true
global IMMUNITY_FROM_DIST = true
global NETWORK_TEST = false 

if IMMUNITY
    global immunity_label = "immunity_on"
else
    global immunity_label = "immunity_off"
end

# 2022 09 19 : tests 
# (2) check that no agents who are removed are tested [x]
# iterate through outbreak-specific input files: 

global data_dir_L1 = pwd() #C:\\Users\\czachreson\\Desktop\\policy_work\\RACF_2022\\RACF_code\\population_generator\\2022_09_02"

global data_dir_L2 = "\\input"#"\\network_generator_step_4\\agents_OB_model_ready_exemplar\\2022_09_02"

global output_dir_L1 = pwd() * "\\test_R0_stratified_contact\\$immunity_label"
if !ispath(output_dir_L1)
    mkpath(output_dir_L1)
end

# NOTE: 2022 09 25 moving outbreak iterator from main() 
# this is the location of the list of outbreak-specific population info
# the population files themselves are contained in the corresponding input directory
global OB_list_dir = "$data_dir_L1\\$data_dir_L2"
global OB_list_fname = "$OB_list_dir\\Outbreak_list_group_exemplars_2022_09_25.csv"
global OB_list = DataFrame(CSV.File(OB_list_fname, delim = ","))
global n_populations = size(OB_list, 1)

# these are the indices of the selected exemplar outbreaks 
# from the initial list of ~60 generated on 2022 08 30
# the indices are retained in the filenames. 
global OB_indices =[6, 28, 57, 33, 55] 

#initialise output file strings, these 
# get re-valued in the scenario loop. 
global delay_label = ""
global t_str = ""
global test_compliance_label = ""

#more globals for iterating through outbreaks
# and generating input/output directory names 
global OB_index = 0 # placeholder (this is the index in the original set of all exemplar outbreaks)
global OB_i = 0 #placeholder (this is the local index for the subset of outbreaks analyzed)
global group_label = ""
global start_date_label = ""
global OB_label = ""
global data_dirname = ""
global output_dir_L2 = ""
global output_dir_OB = ""

# define the main function 
function main()

    
    include("./setup_R0_RACF.jl")

    include("./networks_R0_RACF.jl")

    include("./Diseases_R0_RACF.jl")

    include("./Agents_R0_RACF.jl")

    include("./Facility_Structure_R0.jl")

    include("./Transmission_Dynamics_R0.jl")

    include("./Outbreak_Response_R0.jl")
    

    #output linelist vectors 

    tot_infections = Vector{Int64}()
    tot_detections = Vector{Int64}()
    time_to_detection = Vector{Float64}()
    tot_infections_at_detection = Vector{Float64}()
    generation_interval = Vector{Float64}()

    output_linelist = DataFrame(
        run_id = Int64[], 
        facility_id = Int64[],
        OB_start_date = String[],
        n_residents = Int64[],
        n_staff = Int64[],
        total_FTE = Float64[],
        index_case_type = String[],
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
        FTE_def_max = Float64[]
    )

    #n_outbreaks_tot = 1000 # setting this a global. 
    n_outbreaks = 0
    i = 0
    while i < n_runs
        #i in 1:n_runs
        i += 1
        #println("\n****run:$i")

        global ACTIVE_OUTBREAK = false 

        global seed_offset = i # used for re-seeding rngs for each run, for comparability between each linelist element. 

        apply_seed_offset(seed_offset) #should re-seed random number generators each with a new run-specific value. 
        # this will allow direct comparison of linelist entries from different scenarios on the same facility. 


        # fill any constant params into linelist: 
        push!(output_linelist.run_id, i)
        push!(output_linelist.facility_id, OB_list.service_id[OB_i])
        push!(output_linelist.OB_start_date, OB_list.first_date[OB_i])

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
        
        # 2022 09 14 no longer using this - generating network from room assignments locally 
        #parse network input as neighbour lists:
        #N_lists = N_list()
        #populate_neighbour_lists_from_DataFrame!(N_lists.id_to_contacts, N_lists_str)
        ##
        #initialise Agents 
        agents = Agents()

        if IMMUNITY_FROM_DIST

            ln_neut_dist_res = Normal(log_mu_res, log_sig_res)
            ln_neut_dist_wG = Normal(log_mu_wG, log_sig_wG)
            ln_neut_dist_WM = Normal(log_mu_wM, log_sig_wM)

            populate_workers_from_DataFrame_imDist!(agents.workers_G, workers_G_str, ln_neut_dist_wG)#, N_lists.id_to_contacts)
            N_workers_G = length(agents.workers_G) #note size() not defined for dict
            #println("created $N_workers_G general staff members")

            #initialise medical staff
            populate_workers_from_DataFrame_imDist!(agents.workers_M, workers_M_str, ln_neut_dist_WM)#, N_lists.id_to_contacts)
            N_workers_M = length(agents.workers_M) #note size() not defined for dict
            #println("created $N_workers_M medical staff members")

            #initialise residents
            populate_residents_from_DataFrame_imDist!(agents.residents, residents_str, ln_neut_dist_res)#, N_lists.id_to_contacts)
            N_residents = length(agents.residents) #note size() not defined for dict
            #println("created $N_residents residents")



        else
            # iterate through parsed input DataFrames and initialise workers/residents:
            #initialise general staff
            #NOTE: 2022 09 14 removed initialisation of network from file (done locally from room assignments now)
            populate_workers_from_DataFrame!(agents.workers_G, workers_G_str)#, N_lists.id_to_contacts)
            N_workers_G = length(agents.workers_G) #note size() not defined for dict
            #println("created $N_workers_G general staff members")

            #initialise medical staff
            populate_workers_from_DataFrame!(agents.workers_M, workers_M_str)#, N_lists.id_to_contacts)
            N_workers_M = length(agents.workers_M) #note size() not defined for dict
            #println("created $N_workers_M medical staff members")

            #initialise residents
            populate_residents_from_DataFrame!(agents.residents, residents_str)#, N_lists.id_to_contacts)
            N_residents = length(agents.residents) #note size() not defined for dict
            #println("created $N_residents residents")

        end



        agents.All = merge(agents.residents, agents.workers_G, agents.workers_M)

        # add agent numbers to output list 
        push!(output_linelist.n_residents, N_residents)
        push!(output_linelist.n_staff, (N_workers_G + N_workers_M))
        push!(output_linelist.total_FTE, compute_total_FTE(agents))

        #populate rooms (new on 2022 09 13, need it for re-distribution of labour during 
        # surge rostering from furloughs)
        rooms = Rooms()
        populate_Rooms_from_Agents!(rooms, agents)

        # assign neighbour lists from rooms:
        N_lists_new = N_list_temporal_multigraph() 
        n_days = 7 #hacked in.
        for d in 1:n_days
            N_lists_new.day_to_N_list[d] = N_list()
            N_list_d = N_lists_new.day_to_N_list[d]
            populate_N_lists_d_from_Rooms!(rooms, N_list_d, d)
        end

        # assign agent contacts from N lists constructed above from rooms 
        assign_contacts!(agents, N_lists_new) 
        #note - this will completely replace the existing contact list 

        # adjust pairwise weights based on model-specific bias 
        # NOTE: this involves GLOBAL control parameters. 
        # NOTE: adjusting weights in agent neighbour lists 
        # also adjusts the weights in N_lists (referencing)
        pairwise_weights!(agents)

        if NETWORK_TEST
            # TEST: 2022 09 20
            # copy the network, and test against final network
            # after simulation is complete and all workers 
            # have returned from furlough
            E_list_initial =  E_list()
            for (id, a) in agents.All
                for (d, c) in a.contacts
                    add_source_edges_to_E_list!(E_list_initial, a, c )
                end
            end
        end

        #initialise test schedule 
        initialise_baseline_testing_schedule!(agents)

        ##*****##

        #println("adjusted pairwise weights")

        edge_lists_all = E_list_temporal()
        # iterate through agents and contacts and generate 
        # full edge lists (all possible contacts)
        fill_E_lists_all!(edge_lists_all, agents)
        #println("created edge lists for each day")

        # compute weight totals for each day (denominators used to compute sampling rate)
        w_tot_d = Dict{Int64, Float64}() # day -> w_tot
        for (d, EL) in edge_lists_all.day_to_E_list
            #w_tot_d[d] = sum_edge_weights_EList(EL) # directed
            w_tot_d[d] = sum_edge_weights_EList(EL) / 2.0 #undirected
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
        all_transmissions = DataFrame(source = [], target = [], time = [], type_of_agent = [])

        all_detections = DataFrame(id = [], time_detected = [], type_of_agent = [])
        
        # initialise transmission simulation

        # initialise time keepers 
        t_i = 0.0
        t_f = 90.0 # three months, maximum simulation time. 
        termination_flag = false 
        step_final = convert(Int64, ceil((t_f - t_i)/dt ))

        # TODO: should set the system up to start on an arbitrary day. 
        day_of_week_i = convert(Int64, floor(mod(t_i, 7))) + 1 # first day. 

        #dt is defined in setup # dt = 0.1 (2022 08 17)

        #index case id
        #index_case_id = minimum(keys(agents.workers_G)) #first worker (id 101 in this example) 

        if resident_index_case && worker_index_case
            index_case_id = select_random_agent(agents, day_of_week_i)
        elseif resident_index_case
            index_case_id = select_random_resident(agents)
        elseif worker_index_case
            #index_case_id = select_random_general_staff(agents, day_of_week_i)
            index_case_id = select_random_worker(agents, day_of_week_i)
        end

        if is_worker(agents.All[index_case_id])
            push!(output_linelist.index_case_type, "worker")
        else
            push!(output_linelist.index_case_type, "resident")
        end



        time_infected = t_i

        infect_agent!(agents.All[index_case_id], 
                    diseases["Omicron"], 
                    time_infected, 
                    infected_agents)


        # define some parameters defining contact rates 
        # 2022 09 29 contact rate per high needs resident vs rate per regular-needs resident 
        n_residents_high_needs = count_high_needs_residents(agents)
        n_residents_reg_needs = count_reg_needs_residents(agents)

        contact_rate_per_day = convert(Float64, n_residents_high_needs) * contact_rate_per_resident_per_day
        # applying the same factor of 3 applied to the contact weights between workers and residents. 
        contact_rate_per_day += convert(Float64, n_residents_reg_needs) * (contact_rate_per_resident_per_day / 3.0)

        # background contacts between residents (i.e., in shared meal spaces)
        # these are not included in the structured facility model
        # the implied structure is a fully-connected graph with uniform edge weight
        #bkg_contact_rate_per_resident_per_day = 5.0 #NOTE: moved to setup_RACF.jl
        bkg_contact_rate_per_resident_per_step = bkg_contact_rate_per_resident_per_day * dt

        # 2022 09 23 : including a delay between outbreak declaration and the implementation of measures 
        time_since_outbreak_declared = 0
        global PPE_AVAILABLE = false #this gets toggled on after delay. 

        #for numeric times:
        small_num = 0.0000000001
        step = 0
        while !termination_flag
            step += 1
            
            if step == step_final
                termination_flag = true
            end

            t = convert(Float64, step) * dt 
            day_of_week = convert(Int64, ceil(mod(t-small_num, 7))) # day is an integer 
            day = convert(Int64, ceil(t - small_num)) # first day. 


            #println("time is $t")
            #println("day is $day")

            t_last = convert(Float64, step-1) * dt 
            day_last = convert(Int64, ceil(t_last - small_num)) # day of previous step

            new_day = false
            if day_last < day
                new_day = true
                #println("new day, day of week: $day_of_week")
            end


            # iterate through infected agents and update
            # infection status of each.  
            update_infections!(agents, infected_agents)

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
                update_absentees!(agents, removal_period, removed_workers, t, rooms, room_ids_to_update)
                # room_ids_to_update will now include all rooms where workers have been
                # reinstated from furlough 

                #iterate isolation periods for residents and re-instate
                update_isolated_residents!(agents, removal_period, isolated_residents, t)
                
                F_t = length(removed_workers)
                if F_t > F_max
                    F_max = F_t # maximum number of concurrently furloughed staff 
                end

                Iso_t = length(isolated_residents)
                if Iso_t > Iso_max
                    Iso_max = Iso_t # maximum number of isolated residents 
                end
                
                FTE_def_t = compute_FTE_deficit(agents, removed_workers)
                if FTE_def_t > FTE_def_max
                    FTE_def_max = FTE_def_t # maximum FTE deficit from furlough 
                end
                

                worker_ids_to_remove = Array{Int64, 1}()
                resident_ids_to_isolate = Array{Int64, 1}() 

                # check to see whether to implement infection control after delay: 
                if time_since_outbreak_declared >= delay_val_i
                    global PPE_AVAILABLE = true #toggled on after delay. #see trigger for transmit_infection_AOB in #transmission_dynamics.jl  
                end

                # update delay clock for introduction of infection control measures: 
                if ACTIVE_OUTBREAK
                    time_since_outbreak_declared += 1 # integer increase b/c it's in the new_day scope 
                end

                


                # NOTE: turning outbreak control off for initial R0 test. 
                if (ACTIVE_OUTBREAK && OUTBREAK_CONTROL)
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    test_workers!(all_detections, p_test_per_day_workers_outbreak, agents, infected_agents, day_of_week, t, worker_ids_to_remove)
                    # test residents 
                    test_residents!(all_detections, p_test_per_day_residents_outbreak, agents, infected_agents, day_of_week, t, resident_ids_to_isolate)
                else
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    test_workers!(all_detections, p_test_per_day_workers_baseline, agents, infected_agents, day_of_week, t, worker_ids_to_remove)
                    # test residents 
                    test_residents!(all_detections, p_test_per_day_residents_baseline, agents, infected_agents, day_of_week, t, resident_ids_to_isolate)
                end

                # remove workers from rooms and from N_lists if they tested positive 
                if WORKER_CASE_ISOLATION
                    remove_workers!(agents, worker_ids_to_remove, N_lists_new, removed_workers, t, rooms, room_ids_to_update)
                end

                # isolate residents who tested positive: 
                if RESIDENT_CASE_ISOLATION
                    isolate_residents!(agents, resident_ids_to_isolate, isolated_residents, t)
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
                        update_N_lists_d_from_Rooms!(rooms, room_ids_to_update, N_list_d, d)
                    end
                end
                # assign agent contacts from N lists constructed above from rooms 
                # only doing this for the current day, as the future days will probably change again 
                # before the contacts are needed. 
                assign_todays_contacts!(agents, N_lists_new, day_of_week) 
                #note - this will completely replace the existing contact list for today
        
                # adjust pairwise weights based on model-specific bias 
                pairwise_weights_d!(agents, day_of_week) 
                # for the daily network update, this should be 
                # incorporated into assign_todays_contacts!()

                
                # count detections from last 7 days, if > 0, we have an active outbreak, otherwise we don't_detected
                #println("\n****\nt=$t: new cases last 7 days: $n_cases_last_7")

                # declare active outbreak, or delcare active outbreak over. 
                if !(ACTIVE_OUTBREAK)
                    OB_flag = declare_outbreak(all_detections, t)
                    if OB_flag
                        global ACTIVE_OUTBREAK = true
                        n_outbreaks += 1 
                        #println("active outbreak declared at time $t")
                        #println("there are: $(size(all_transmissions,1)) infections ")
                        output_linelist.OB_declared[i] = true 
                        output_linelist.t_OB_on[i] = t
                        output_linelist.I_cum_OB_on[i] = size(all_transmissions,1) + 1 # add one to include index case 
                        output_linelist.Det_cum_OB_on[i] = size(all_detections, 1) 
                    end
                else
                    OB_over_flag = outbreak_over(all_detections, t)
                    if OB_over_flag
                        global ACTIVE_OUTBREAK = false
                        #println("active outbreak over at time $t")
                        output_linelist.t_OB_off[i] = t
                        termination_flag = true
                    end
                end

            end 
            
            # define some parameters defining contact rates 
            #contact_rate_per_resident_per_day = 10.0 #NOTE: moved to setup_RACF.jl

            # 2022 09 28 contact rate per high needs resident vs rate per regular-needs resident 
            n_residents_high_needs = count_high_needs_residents(agents)
            n_residents_reg_needs = count_reg_needs_residents(agents)

            contact_rate_per_day = convert(Float64, n_residents_high_needs) * contact_rate_per_resident_per_day
            # applying the same factor of 3 applied to the contact weights between workers and residents. 
            contact_rate_per_day += convert(Float64, n_residents_reg_needs) * (contact_rate_per_resident_per_day / 3.0)

            #contact_rate_per_day = convert(Float64, N_residents) * contact_rate_per_resident_per_day
            
            contact_rate_per_step = contact_rate_per_day * dt 

            # background contacts between residents (i.e., in shared meal spaces)
            # these are not included in the structured facility model
            # the implied structure is a fully-connected graph with uniform edge weight
            # bkg_contact_rate_per_resident_per_day = 5.0 #NOTE: moved to setup_RACF.jl

            bkg_contact_rate = bkg_contact_rate_per_resident_per_step
            bkg_contact_rate_iso = bkg_contact_rate_per_resident_per_step * (1.0 - resident_isolation_efficacy)

            if ACTIVE_OUTBREAK && OUTBREAK_CONTROL && RESIDENT_LOCKDOWN
                bkg_contact_rate *= (1.0 - resident_lockdown_efficacy)
            end


            compute_transmission_R0!(index_case_id, all_transmissions, infected_agents, agents,
                                day_of_week, w_tot_d, contact_rate_per_step, bkg_contact_rate, bkg_contact_rate_iso, t)

        end

        # 2022 09 20
        #TEST: build final network and compare to original one
        # (should be the same if all workers have been reinstated by tf)
        if NETWORK_TEST
            E_list_final =  E_list()
            for (id, a) in agents.All
                for (d, c) in a.contacts
                    add_source_edges_to_E_list!(E_list_final, a, c )
                end
            end

            test_val_e1_e2 = compare_E_lists(E_list_initial, E_list_final)
            if test_val_e1_e2 
                println("final network same as initial network e_i -> e_f")
            else
                println("final network different from initial network e_i -> e_f")
            end
            test_val_e2_e1 = compare_E_lists(E_list_final, E_list_initial)
            if test_val_e2_e1 
                println("final network same as initial network e_f -> e_i")
            else
                println("final network different from initial network e_f -> e_i")
            end
        end

        if write_transmission_tree_flag
            rand_label = now()
            rand_label = replace(string(rand_label), ":"=> "_")
            rand_label = replace(string(rand_label), "."=> "_")
            rand_label = replace(string(rand_label), "-"=> "_")
            if OUTBREAK_CONTROL
                ob_lab = "ICAOB_"
            else
                ob_lab = "unmitigated_"
            end
            
            output_dir_run = "$output_dir_OB\\run_$(i)"
            if !ispath(output_dir_run)
                mkpath(output_dir_run)
            end

            CSV.write("$(output_dir_run)\\$(ob_lab)detected_cases_$rand_label.csv", all_detections)
            CSV.write("$(output_dir_run)\\$(ob_lab)transmission_tree_$rand_label.csv", all_transmissions)
        end

        push!(output_linelist.F_max, F_max)# maximum number of concurrently furloughed staff 
        push!(output_linelist.Iso_max, Iso_max) # maximum number of isolated residents 
        push!(output_linelist.FTE_def_max, FTE_def_max) # maximum FTE deficit from furlough 


        # count total infections and detections in staff and residents 

        I_tot_staff = count_worker_infections(agents, all_transmissions)
        I_tot_res = count_resident_infections(agents, all_transmissions)
        Det_tot_staff = count_worker_detections(agents, all_detections)
        Det_tot_res = count_resident_detections(agents, all_detections)

        if is_worker(agents.All[index_case_id])
            I_tot_staff += 1
        else 
            I_tot_res += 1
        end

        push!(output_linelist.I_tot_staff, I_tot_staff) 
        push!(output_linelist.Det_tot_staff, Det_tot_staff) 
        push!(output_linelist.I_tot_res, I_tot_res) 
        push!(output_linelist.Det_tot_res, Det_tot_res) 

        t_gen = NaN 



        if isempty(all_transmissions)
            total_infections = 1
            
        else
            total_infections = size(all_transmissions, 1) + 1 # + 1 is for index case. 
            t_gen = mean(all_transmissions.time)
        end

        #println("T_gen for run: $t_gen")
        push!(generation_interval,t_gen)

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
        push!(time_to_detection, t_first_detection)
        push!(tot_detections, total_detections)
        push!(tot_infections_at_detection, total_infections_at_detection)
    
    end

    println("avg. total infections: $(mean(tot_infections))")
    println("avg. total detected cases: $(mean(tot_detections))")
    println("avg. t first detection: $(mean(time_to_detection[.!isnan.(time_to_detection)]))")
    println("avg. infections at first detection: $(mean(tot_infections_at_detection[.!isnan.(tot_infections_at_detection)]))")
    println("avg. generation interval:  $(mean(generation_interval[.!isnan.(generation_interval)]))")

    output_fname = "$(output_dir_OB)\\output_summary.csv"
    output_linelist_fname = "$(output_dir_OB)\\output_linelist.csv"

    output = DataFrame(I_tot = tot_infections, 
                       I_detected = tot_detections, 
                        t_first_detection = time_to_detection, 
                        I_t_detection = tot_infections_at_detection,
                        T_gen = generation_interval)

    CSV.write(output_fname, output)
    CSV.write(output_linelist_fname, output_linelist)



end



# loop through exemplars and run the simulations: 

for i in 1:n_populations#outbreak_indices

    global OB_i = i #OB_i gets used in main() for setting up outputs. 
    global OB_index = OB_indices[OB_i]
    global group_label = "group_$(OB_list.group_id[OB_i])"
    global start_date_label = replace(OB_list.first_date[OB_i], "/"=> "_")
    global OB_label = "OB_$(OB_index)_facID_$(OB_list.service_id[OB_i])_$start_date_label"
    global data_dirname = "$(data_dir_L1)$(data_dir_L2)\\$group_label\\$OB_label"#"$(data_dir_L1)$(data_dir_L2)\\$group_label\\$OB_label"

    if !ispath(data_dirname)
        continue 
    end


    for (d_i) in delay_vals 
        global delay_val_i = d_i
        for (t_j) in test_compliance
            global test_compliance_j = t_j


            global delay_label = "response_delay_$delay_val_i"
            global t_str = replace("$test_compliance_j", "."=> "p")
            global test_compliance_label = "test_compliance_$t_str"

            global output_dir_L2 = "$group_label\\$OB_label\\$delay_label\\$test_compliance_label"
            global output_dir_OB = "$(output_dir_L1)\\$(output_dir_L2)"
            if !ispath(output_dir_OB)
                mkpath(output_dir_OB)
            end

            main()

        end
    end
end