
module Outbreak_Response

using StatsBase
using Random
using DataFrames

using Main.Setup_RACF
using Main.Agents_RACF
using Main.Networks_RACF
using Main.Facility_Structure
using Main.Diseases_RACF


# simulate testing and detection 

function initialise_baseline_testing_schedule!(agents::Agents_RACF.Agents_T, config::Setup_RACF.Config_T) #agents 
    # no baseline testing schedule for residents 
    # try to test workers once per roster cycle if part-time, and twice per roster schedule if full time 
    # (days 1 and 4)
    for (id, a) in agents.workers_G 
        assign_test_days!(a, config)
    end

    for (id, a) in agents.workers_M
        assign_test_days!(a, config)
    end

end

function assign_test_days!(a::Agents_RACF.Agent_T, config::Setup_RACF.Config_T)
    r = a.roster
    a.test_schedule = zeros(Int64, size(a.roster))
    PT_flag = false
    # part time?
    if sum(r) < 5
        PT_flag = true 
    else
        PT_flag = false 
    end

    if PT_flag 
    # locate work days
    days_working = findall(==(1), a.roster)
    test_day = sample(config.rng_testing, days_working) # config
    a.test_schedule[test_day] = 1
    else
        days_working = findall(==(1), a.roster)
        

        test_day_1 = sample(config.rng_testing, days_working) # config
        test_day_2 = mod1(test_day_1 + 3, 7)

        while (a.roster[test_day_2] == 0) || (abs(test_day_1 - test_day_2) < 3)
            test_day_1 = sample(config.rng_testing, days_working) # config
            test_day_2 = mod1(test_day_1 + 3, 7)
        end
        
        a.test_schedule[test_day_1] = 1
        a.test_schedule[test_day_2] = 1
    end

end

function test_agent(a::Agents_RACF.Agent_T, config::Setup_RACF.Config_T)::Bool

    positive = false 

    for (pathogen_name, infection) in a.infections
        Diseases_RACF.compute_test_sensitivity!(infection)
        if rand(config.rng_testing) < infection.test_sensitivity #config
            positive = true
        end
    end

    return positive
end

function test_workers!(detections::DataFrame, p_test_per_day::Float64, 
                       agents::Agents_RACF.Agents_T, infected_agents::Dict{Int64, Float64},
                       day_of_week::Int64, t::Float64, ids_to_remove::Array{Int64, 1},
                       config::Setup_RACF.Config_T) # agents, df, 


    # test for existing infections (RAT): 
    for (id, t_inf) in infected_agents
        # check if they've already been detected: 
        
        a = agents.All[id] 
        if a.t_detected < 0.0
            if Agents_RACF.is_worker(a)
            
                test_this_worker = false 
                if !config.active_outbreak
                    test_this_worker = ( (a.test_schedule[day_of_week] == 1) && (a.t_removed < 0.0 ) )
                else 
                    test_this_worker = ((a.roster[day_of_week] == 1) && (a.t_removed < 0.0 ) )# if they're present 
                end

                if test_this_worker# if they're present and scheduled for a test 
                    
                    p = p_test_per_day
                    if check_symptom_expression(a)
                        p = config.p_test_if_symptomatic
                    end
                    
                    if rand(config.rng_testing) < p #config
                        tested_positive = test_agent(a, config)
                        if tested_positive 

                            #println("Agent $(a.id) test schedule: $(a.test_schedule), day of week: $day_of_week")

                            a.t_detected = t

                            push!(ids_to_remove, a.id)
                            
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
                         agents::Agents_RACF.Agents_T, infected_agents::Dict{Int64, Float64}, 
                         day_of_week::Int64, t::Float64, resident_ids_to_isolate::Array{Int64, 1},
                         config::Setup_RACF.Config_T) #agents

    # test for existing infections (RAT): 
    for (id, t_inf) in infected_agents
        a = agents.All[id] 
        # NOTE: if we change the model to SIS, this will need to be modified. 
        # the variable t_removed gets reinitialisd to -1.0 after iso completes. 
        if a.t_detected < 0.0 # initialised negative for this purpose 
            if Agents_RACF.is_resident(a)
            
                if haskey(a.contacts, day_of_week) # if they're present 

                    p = p_test_per_day
                    if check_symptom_expression(a)
                        p = config.p_test_if_symptomatic
                    end
                    
                    if rand(config.rng_testing) < p #config
                        tested_positive = test_agent(a, config)
                        if tested_positive 

                            push!(resident_ids_to_isolate, a.id)

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


function check_symptom_expression(a::Agents_RACF.Agent_T)::Bool #agents

    symptoms = false 
    for (name, im) in a.infections
        if im.expressing_symptoms
            symptoms = true
        end
    end

    return symptoms

end

# toggle active outbreak

# active outbreak definition: 
 # 2 or more resident cases in 5 days 
 # OR 
 # 5 or more worker cases in 7 days 

function count_detections_last_7_days(detections::DataFrame, t)::Int64
   
    d = floor(t) 
    dm7 = d - 7

    detections_last_7 = detections[detections.time_detected .> dm7, :]

    return size(detections_last_7, 1)

end

function count_resident_detections_last_5_days(detections::DataFrame, t)::Int64
   
    d = floor(t) 
    dm5 = d - 5

    #res_detections_last_5 = subset(detections, :time_detected => ByRow(>(dm5)), :type_of_agent => ByRow(==(1)))
    res_detections_last_5 = detections[((detections.time_detected .> dm5) .& (detections.type_of_agent .== 1)), :]
    # TODO: should make agent typing system in detections list more formal. 
    # (currently  1 = resident, 2 = gen staff, 3 = med staff)

    return size(res_detections_last_5, 1)

end

function count_resident_detections_last_7_days(detections::DataFrame, t)::Int64
   
    d = floor(t) 
    dm7 = d - 7

    res_detections_last_7 = detections[((detections.time_detected .> dm7) .& (detections.type_of_agent .== 1)), :]
    # TODO: should make agent typing system in detections list more formal. 
    # (currently  1 = resident, 2 = gen staff, 3 = med staff)

    return size(res_detections_last_7, 1)

end

function count_worker_detections_last_7_days(detections::DataFrame, t)::Int64
   
    d = floor(t) 
    dm7 = d - 7

    staff_detections_last_7 = detections[((detections.time_detected .> dm7) .& (detections.type_of_agent .> 1)), :]

    return size(staff_detections_last_7, 1)

end

function declare_outbreak(detections::DataFrame, t)::Bool

    OB_on = false

    res_last_5 = count_resident_detections_last_5_days(detections, t)
    staff_last_7 = count_worker_detections_last_7_days(detections, t) 

    if res_last_5 > 2 || staff_last_7 > 5
        OB_on = true 
    end

    return OB_on 
end

# declaring an outbreak over, 7 days of no new resident cases 
function outbreak_over(detections::DataFrame, t)::Bool

    OB_over = false 
    res_cases_7 = count_resident_detections_last_7_days(detections, t)
    if res_cases_7 == 0
        OB_over = true 
    end

    return OB_over 
end

# remove a worker : NOTE: don't need to modify w_tot_d IF we assume that removing a worker
# reduces the overall contact rate proportionally to w_i / w_tot (differences cancel) [w_i is weight of agent i's contacts]

#RE-DISTRIBUTION 
# when a worker is removed, some of their contact weight gets distributed 
# to other long-term employees present, and some gets distributed to a 
# temporary employee (these portions can be altered as control parameters)

function remove_worker!(a::Agents_RACF.Agent_T, 
                        t_removed::Float64, 
                        removed_agents::Dict{Int64, Float64}) #agents
    # a removed worker can't infect or be infected. 
    a.t_removed = t_removed
    removed_agents[a.id] = t_removed

end

function remove_worker_from_rooms!(a::Agents_RACF.Agent_T, 
                                   rooms::Facility_Structure.Rooms_T, 
                                   room_ids_to_update::Set{Int64}) #agents, facility
    
    # iterate through agent's rooms and remove them from the appropriate sets: 
    for (d, room_list) in a.rooms
        for rm_id in room_list
            room = rooms.Day_to_Rooms[d][rm_id]
            delete!(room.agent_ids, a.id )
            delete!(room.worker_ids, a.id)

            if !in(rm_id, room_ids_to_update)
                push!(room_ids_to_update, rm_id)
            end

            room.N_workers_t -= 1
        end
    end
 
end

function remove_workers!(agents::Agents_RACF.Agents_T, 
                         ids_to_remove::Array{Int64, 1}, 
                         N_list_TMG::Networks_RACF.N_list_temporal_multigraph_T, 
                         removed_workers::Dict{Int64, Float64}, 
                         t::Float64, 
                         rooms::Facility_Structure.Rooms_T, 
                         room_ids_to_update::Set{Int64}) #agents, networks, facility

    #remove worker agent: 
    for id in ids_to_remove
        a = agents.All[id]
        remove_worker!(a, t, removed_workers)
        remove_worker_from_rooms!(a, rooms, room_ids_to_update)

        for (d, N_list_d) in N_list_TMG.day_to_N_list
        N_list_d.id_to_contacts[a.id] = Array{Networks_RACF.Contact_T, 1}()  #networks 
        # clear all contacts associated with removed worker (out edges). #
        # NOTE: in-edges will be removed from contact's edge lists during the room-by-room network update.  
        end

    end


end


function isolate_resident!(a::Agents_RACF.Agent_T, 
                           t_isolated::Float64, 
                           isolated_residents::Dict{Int64, Float64}) #agents
    # isolated residents still interact with staff and roommates, but do not participate in random interactions
    # with other residents 
    a.t_removed = t_isolated
    isolated_residents[a.id] = t_isolated
end

function isolate_residents!(agents::Agents_RACF.Agents_T, 
                            ids_to_isolate::Array{Int64, 1}, 
                            isolated_residents::Dict{Int64, Float64}, 
                            t::Float64) #agents 

    for id in ids_to_isolate
        a = agents.residents[id]
        #isolate resident (if RESIDENT_CASE_ISOLATION flag is true)
        isolate_resident!(a, t, isolated_residents)
    end

end



function update_absentees!(agents::Agents_RACF.Agents_T, 
                           removal_period::Float64, 
                           removed_agents::Dict{Int64, Float64}, 
                           t::Float64, 
                           rooms::Facility_Structure.Rooms_T, 
                           room_ids_to_update::Set{Int64} ) #agents, facility

    n_absent = 0 

    agents_to_reinstate = []
    for (id_i, t_removed) in removed_agents
        if t - t_removed > removal_period
            a = agents.All[id_i]
            a.t_removed = -1.0 #not removed 
            push!(agents_to_reinstate, id_i)
        else
            n_absent += 1
        end
    end

    # add back to rooms (I know, I know, assignment of a is redundant...)
    for id_i in agents_to_reinstate

        a = agents.All[id_i]

        # 2022 09 15: adding re-instated agents back to rooms. 
        
        # iterate through agent's rooms and add them to the appropriate sets: 
        for (d, room_list) in a.rooms
            for rm_id in room_list
                room = rooms.Day_to_Rooms[d][rm_id]
                push!(room.agent_ids, a.id )
                push!(room.worker_ids, a.id)
                if !in(rm_id, room_ids_to_update)
                    push!(room_ids_to_update, rm_id)
                end
                room.N_workers_t -= 1
            end
        end 

        delete!(removed_agents, id_i)
    end

    #println("number of agents absent: $n_absent")

end
# reinstate a worker (modifies the contact weight denominators)

function update_isolated_residents!(agents::Agents_RACF.Agents_T, 
                                    removal_period::Float64, 
                                    isolated_residents::Dict{Int64, Float64}, 
                                    t::Float64) #agents 

    n_isolated = 0 

    agents_to_release = []
    for (id_i, t_removed) in isolated_residents
        if t - t_removed > removal_period
            a = agents.All[id_i]
            a.t_removed = -1.0 #not isolated 
            push!(agents_to_release, id_i)
        else
            n_isolated += 1
        end
    end

    for id_i in agents_to_release
        delete!(isolated_residents, id_i)
    end

    #println("number of agents absent: $n_isolated")

end

function is_isolated(a::Agents_RACF.Agent_T)::Bool #agents 
    iso_flag = false 
    if a.t_removed >= 0.0
        iso_flag  = true
    end
    return iso_flag 
end




end