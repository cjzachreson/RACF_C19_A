
#TODO: update_infection!() should compute FoI, symptom expression, test sensitivity, and recovery. 

include("./header_RACF.jl")

# global control parameter: worker index case vs. resident index case (if both are true it can be either)
global resident_index_case = false
global worker_index_case = true

if resident_index_case && worker_index_case
    index_case_label = "resident_or_worker_index"
elseif resident_index_case
    index_case_label = "resident_index"
elseif worker_index_case
    index_case_label = "worker_index"
end

# iterate through outbreak-specific input files: 

data_dir_L1 = "C:\\Users\\czachreson\\Desktop\\policy_work\\RACF_2022\\RACF_code\\population_generator\\2022_09_02"

data_dir_L2 = "\\network_generator_step_4\\agents_OB_model_ready_exemplar\\2022_09_02"

output_dir_L1 = pwd() * "\\exemplar_outbreak_results"
if !ispath(output_dir_L1)
    mkpath(output_dir_L1)
end

OB_list_dir = "$data_dir_L1\\Outbreak_list_step_1"
OB_list_fname = "$OB_list_dir\\Outbreak_list_g2_g5_exemplars_2022_09_05.csv"
OB_list = DataFrame(CSV.File(OB_list_fname, delim = ","))

n_outbreaks = size(OB_list, 1)

outbreak_indices =[17, 43] #HACK 

println("$n_outbreaks")

for OB_i in 1:n_outbreaks#outbreak_indices

    outbreak_index = outbreak_indices[OB_i]
    
    group_label = "group_$(OB_list.group_id[OB_i])"

    start_date_label = replace(OB_list.first_date[OB_i], "/"=> "_")

    OB_label = "OB_$(outbreak_index)_facID_$(OB_list.service_id[OB_i])_$start_date_label"

    global data_dirname = "$(data_dir_L1)$(data_dir_L2)\\$group_label\\$OB_label"

    if !ispath(data_dirname)
        continue 
    end


    output_dir_L2 = "$index_case_label\\$group_label\\$OB_label"

    global output_dir_OB = "$(output_dir_L1)\\$(output_dir_L2)"

    if !ispath(output_dir_OB)
        mkpath(output_dir_OB)
    end


    include("./setup_RACF.jl")

    include("./networks_RACF.jl")

    include("./Diseases_RACF.jl")

    include("./Agents_RACF.jl")

    include("./Transmission_Dynamics.jl")

    tot_infections = Vector{Int64}()
    tot_detections = Vector{Int64}()
    time_to_detection = Vector{Float64}()
    tot_infections_at_detection = Vector{Float64}()


    #@time begin
    n_runs = 1000

    for i in 1:n_runs

        #println("$i")


        #flag telling the system to write the transmission tree and detections to file: 
        write_output_flag = false

        ## ***** ##
        #parse network input as neighbour lists:

        N_lists = N_list_temporal()
        populate_neighbour_lists_from_DataFrame!(N_lists.id_to_contacts, N_lists_str)
        ##
        #initialise Agents 
        agents = Agents()

        # iterate through parsed input DataFrames and initialise workers/residents:
        #initialise general staff
        populate_workers_from_DataFrame!(agents.workers_G, workers_G_str, N_lists.id_to_contacts)
        N_workers_G = length(agents.workers_G) #note size() not defined for dict
        #println("created $N_workers_G general staff members")

        #initialise medical staff
        populate_workers_from_DataFrame!(agents.workers_M, workers_M_str, N_lists.id_to_contacts)
        N_workers_M = length(agents.workers_M) #note size() not defined for dict
        #println("created $N_workers_M medical staff members")

        #initialise residents
        populate_residents_from_DataFrame!(agents.residents, residents_str, N_lists.id_to_contacts)
        N_residents = length(agents.residents) #note size() not defined for dict
        #println("created $N_residents residents")

        agents.All = merge(agents.residents, agents.workers_G, agents.workers_M)


        ##*****##

        # adjust pairwise weights based on model-specific bias 
        # NOTE: this involves GLOBAL control parameters. 
        # NOTE: adjusting weights in agent neighbour lists 
        # also adjusts the weights in N_lists (referencing)
        pairwise_weights!(agents)
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
        removed_agents = Dict{Int64, Float64}()

        # output structures 
        all_transmissions = DataFrame(source = [], target = [], time = [], type_of_agent = [])

        all_detections = DataFrame(id = [], time_detected = [], type_of_agent = [])
        
        # initialise transmission simulation

        # initialise time keepers 
        t_i = 0.0
        t_f = 90.0 # three months
        step_final = convert(Int64, ceil((t_f - t_i)/dt ))

        day_of_week_i = convert(Int64, floor(mod(t_i, 7))) + 1 # first day. 
        day_i = convert(Int64, floor(mod(t_i, 1))) + 1 # first day. 

        #dt is defined in setup # dt = 0.1 (2022 08 17)

        #index case id
        #index_case_id = minimum(keys(agents.workers_G)) #first worker (id 101 in this example) 

        if resident_index_case && worker_index_case
            index_case_id = select_random_agent(agents)
        elseif resident_index_case
            index_case_id = select_random_resident(agents)
        elseif worker_index_case
            #index_case_id = select_random_general_staff(agents, day_of_week_i)
            index_case_id = select_random_worker(agents, day_of_week_i)
        end



        time_infected = t_i

        infect_agent!(agents.All[index_case_id], 
                      diseases["Omicron"], 
                      time_infected, 
                      infected_agents)


        # define some global parameters defining contact rates 
        contact_rate_per_resident_per_day = 10.0
        contact_rate_per_day = convert(Float64, N_residents) * contact_rate_per_resident_per_day
        contact_rate_per_step = contact_rate_per_day * dt 

        # background contacts between residents (i.e., in shared meal spaces)
        # these are not included in the structured facility model
        # the implied structure is a fully-connected graph with uniform edge weight
        bkg_contact_rate_per_resident_per_day = 5.0
        bkg_contact_rate_per_resident_per_step = bkg_contact_rate_per_resident_per_day * dt


        #for numeric times:
        small_num = 0.0000000001

        for step = 1:step_final
            t = convert(Float64, step) * dt 
            day_of_week = convert(Int64, ceil(mod(t+small_num, 7))) # day is an integer 
            day = convert(Int64, ceil(t)) # first day. 

            t_last = convert(Float64, step-1) * dt 
            day_last = convert(Int64, ceil(t_last)) # day of previous step


            new_day = false
            if day_last < day
                new_day = true
            end


            # iterate through infected agents and update
            # infection status of each.  
            update_infections!(agents, infected_agents)

            #check for symptoms: 

            # test with RATs 
            if new_day 

                #println("\n*****")
                #println("time: $t, day of week: $day_of_week")
                #println( "today is $day, last day was: $day_last")

                update_absentees!(agents, removal_period, removed_agents, t)

                #test_agents!(all_detections, p_test_per_day, agents, infected_agents, day_of_week, t)
                

                if (ACTIVE_OUTBREAK && OUTBREAK_CONTROL)
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    test_workers!(all_detections, p_test_per_day_workers_outbreak, agents, infected_agents, day_of_week, t, removed_agents)
                    # test residents 
                    test_residents!(all_detections, p_test_per_day_residents_outbreak, agents, infected_agents, day_of_week, t)
                else
                    # test workers, workers who test positive are removed for (e.g.) 14 days 
                    test_workers!(all_detections, p_test_per_day_workers_baseline, agents, infected_agents, day_of_week, t, removed_agents)
                    # test residents 
                    test_residents!(all_detections, p_test_per_day_residents_baseline, agents, infected_agents, day_of_week, t)
                end


                # count detections from last 7 days, if > 0, we have an active outbreak, otherwise we don't_detected
                n_cases_last_7 = count_detections_last_7_days(all_detections, t)

                #println("\n****\nt=$t: new cases last 7 days: $n_cases_last_7")

                if n_cases_last_7 > 0
                    global ACTIVE_OUTBREAK = true
                    #println("outbreak active")
                    
                else
                    global ACTIVE_OUTBREAK = false
                end

            end 
            
            # this now implements PPE when outbreak is active. 
            #modified intputs: all_transmissions, infected_agents, agents

            bkg_contact_rate = bkg_contact_rate_per_resident_per_step
            if ACTIVE_OUTBREAK && OUTBREAK_CONTROL && RESIDENT_LOCKDOWN
                bkg_contact_rate *= (1.0 - resident_lockdown_efficacy)
            end

            compute_transmission!(all_transmissions, infected_agents, agents,
                                  day_of_week, w_tot_d, contact_rate_per_step, bkg_contact_rate, t)

        end

        if write_output_flag
            rand_label = now()
            rand_label = replace(string(rand_label), ":"=> "_")
            rand_label = replace(string(rand_label), "."=> "_")
            rand_label = replace(string(rand_label), "-"=> "_")
            if OUTBREAK_CONTROL
                ob_lab = "ICAOB_"
            else
                ob_lab = "unmitigated_"
            end
            
            CSV.write("$(ob_lab)detected_cases_$rand_label.csv", all_detections)
            CSV.write("$(ob_lab)transmission_tree_$rand_label.csv", all_transmissions)
        end

        if isempty(all_transmissions)
            total_infections = 1
        else
            total_infections = size(all_transmissions, 1) + 1 # + 1 is for index case. 
        end
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

    println("total infections: $(mean(tot_infections))")
    println("total detected cases: $(mean(tot_detections))")
    println("t first detection: $(mean(time_to_detection[.!isnan.(time_to_detection)]))")
    println("infections at first detection: $(mean(tot_infections_at_detection[.!isnan.(tot_infections_at_detection)]))")

    output_fname = "$(output_dir_OB)\\output_summary.csv"

    output = DataFrame(I_tot = tot_infections, 
                       I_detected = tot_detections, 
                       t_first_detection = time_to_detection, 
                       I_t_detection = tot_infections_at_detection)

    CSV.write(output_fname, output)
    #TODO: write config. 

end