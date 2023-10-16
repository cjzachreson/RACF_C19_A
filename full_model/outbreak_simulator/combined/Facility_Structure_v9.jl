#Facility_structure.jl

module Facility_Structure

# define room structure: 

using Main.Agents_RACF
using Main.Networks_RACF

abstract type Room_T end 
abstract type Rooms_T end 

abstract type Facility_T end


mutable struct facility<:Facility_T
    id::Int64
    facility() = new()
end


mutable struct Room <: Room_T

    # set or array? don't care about order, but would like to add and
    # delete objects efficiently.
    # NOTE: sets are not ordered, but can be iterated - removing elements by value is fast  
    id::Int64 #unique identifier 
    agent_ids::Set{Int64} 
    worker_ids::Set{Int64}
    resident_ids::Set{Int64}
    N_workers_t0::Int64 # day -> n workers at time t0
    N_workers_t::Int64 # day -> n workers at time t
    N_residents::Int64 # number of resident occupants 

    Room(id) = new( # constructor requires a unique id. 
        id,
        Set{Int64}(), #agent_ids 
        Set{Int64}(), #worker_ids 
        Set{Int64}(), #resident_ids
        0, #N_workers_t0
        0, #N_workers_t
        0  #N_residents 
    )
end


mutable struct Rooms <: Rooms_T 
    # maps day -> [room_id -> room], each room has a set of agents 
    Day_to_Rooms::Dict{Int64, Dict{Int64, Room_T}}
    
    Rooms() = new(
        Dict{Int64, Dict{Int64, Room_T}}()
    )
end




# Initialise rooms sturcture: 

function populate_Rooms_from_Agents!(rooms::Rooms_T, agents::Agents_RACF.Agents_T)


    #identify number of days in roster 
    n_roster_days = 0

    # iterate through agents: 
    for (id, a) in agents.workers_G
        if n_roster_days == 0
            n_roster_days = size(a.roster, 1)
        end
        # for a worker, the room entry is 
        # dict day => set of rooms 
        for (d, room_ids) in a.rooms 
            # check if Rooms has any entries for d 
            if !haskey(rooms.Day_to_Rooms, d)
                # initialise empty room set for day d 
                rooms.Day_to_Rooms[d] = Dict{Int64, Room_T}()
            end
            rooms_d = rooms.Day_to_Rooms[d]
            #now iterate through room ids 
            for rm_id in room_ids 
                # check if room list for day d has room id
                if !haskey(rooms_d, rm_id)
                    rooms_d[rm_id] = Room(rm_id)
                end
                room = rooms_d[rm_id]
                # add worker agent to room 
                push!(room.agent_ids, a.id) 
                push!(room.worker_ids, a.id) 
                room.N_workers_t0 += 1 # day -> n workers at time t0
                room.N_workers_t += 1 # day -> n workers at time t
            end
        end
    end

    for (id, a) in agents.workers_M
        # for a worker, the room entry is 
        # dict day => set of rooms 
        for (d, room_ids) in a.rooms 
            # check if Rooms has any entries for d 
            if !haskey(rooms.Day_to_Rooms, d)
                # initialise empty room set for day d 
                rooms.Day_to_Rooms[d] = Dict{Int64, Room_T}()
            end
            rooms_d = rooms.Day_to_Rooms[d]
            #now iterate through room ids 
            for rm_id in room_ids 
                # check if room list for day d has room id
                if !haskey(rooms_d, rm_id)
                    rooms_d[rm_id] = Room(rm_id)
                end
                room = rooms_d[rm_id]
                # add worker agent to room 
                push!(room.agent_ids, a.id) 
                push!(room.worker_ids, a.id) 
                room.N_workers_t0 += 1 # day -> n workers at time t0
                room.N_workers_t += 1 # day -> n workers at time t
            end
        end
    end

    for (id, a) in agents.residents
        # for a resident, the room is just an integer 
        rm_id = a.room 
        for d in 1:n_roster_days
            # check if Rooms has any entries for d 
            if !haskey(rooms.Day_to_Rooms, d)
                # check room set for day d 
                #Rooms.Day_to_Rooms[d] = Dict{Int64, Room_T}()
                println("WARNING: room list for day $d initialised with no workers")
            end
            rooms_d = rooms.Day_to_Rooms[d]
            
            # check if room list for day d has room id
            if !haskey(rooms_d, rm_id)
                println("WARNING: room $rm_id on day $d initialised with no workers")
            end
            room = rooms_d[rm_id]
            # add worker agent to room 
            push!(room.agent_ids, a.id) 
            push!(room.resident_ids, a.id) 
            room.N_residents += 1
        
        end
    end


end

##### used in R0 implementation
function populate_Rooms_from_Agents_no_workers!(rooms::Rooms_T, agents::Agents_RACF.Agents_T)


    #identify number of days in roster 
    n_roster_days = 7

    # iterate through agents: 

    for (id, a) in agents.residents
        # for a resident, the room is just an integer 
        rm_id = a.room 
        for d in 1:n_roster_days
            # check if Rooms has any entries for d 
            if !haskey(rooms.Day_to_Rooms, d)
                # check room set for day d 
                rooms.Day_to_Rooms[d] = Dict{Int64, Room_T}()
                
            end
            rooms_d = rooms.Day_to_Rooms[d]
            
            # check if room list for day d has room id
            if !haskey(rooms_d, rm_id)
                rooms_d[rm_id] = Room(rm_id)
            end
            room = rooms_d[rm_id]
            # add resident agent to room 
            push!(room.agent_ids, a.id) 
            push!(room.resident_ids, a.id) 
            room.N_residents += 1
        
        end
    end


end


#using N_lists structure from networks_RACF.jl

# modifies N_lists_in based on current room allocations 
# can then call 'add_neighbours_from_N_list(agent, N_lists, agent_id) 
# for each affected agent 
function populate_N_lists_d_from_Rooms!(rooms::Rooms_T, N_lists_in::Networks_RACF.N_list_T, d::Int64)

    # go back to multigraph system... 
    # currently N_lists is initialised as a weighted network, 
    # but this will be harder to update when weights must be aggregated
    # from contacts involving different rooms,
    # add room id to contact object to make multigraph more 
    # explicit  

    rooms_d = rooms.Day_to_Rooms[d]

    id_to_contacts = N_lists_in.id_to_contacts

    #iterate through rooms 
    for (room_id, room) in rooms_d 
        # iterate through agents in room 

        #debugging: 
        #if (room_id == 79 ) && (d == 4)
         #   println("check agent 174")
        #end

        for source_id in room.agent_ids
        # assign out edges to all other agents 
        #(assign each a weight of 1 b/c it's a multigraph)
            for target_id in room.agent_ids
                if target_id != source_id
                    weight = 1.0
                    contact = Networks_RACF.Contact(target_id, d, weight)
                    contact.room_id = room_id 
                    if haskey(id_to_contacts, source_id)
                        push!(id_to_contacts[source_id], contact) 
                        # note: when adapting this to update existing edge lists, I'll have to first clear all old contacts FROM this room_id
                        # I think I can use filter!() and: 
                        # contacts_a = agents.All[source_id].contacts[day_of_week]
                        # filter!((c->c.room_id).(contacts_a) .!= room_id, contacts_a) # this will modify the contacts to remove all those without the current room_id
                    else
                        id_to_contacts[source_id] = Array{Networks_RACF.Contact_T, 1}()
                        push!(id_to_contacts[source_id], contact)
                    end 

                end
            end
        end
    end
end


function update_N_lists_d_from_Rooms!(rooms::Rooms_T, 
                                      room_ids_to_update::Set{Int64}, 
                                      N_lists_in::Networks_RACF.N_list_T, 
                                      d::Int64)
    #(Rooms, room_ids_to_update, N_list_d, day_of_week)

    # should go back to multigraph system... 
    # currently N_lists is initialised as a weighted network, 
    # but this will be harder to update when weights must be aggregated
    # from contacts involving different rooms,
    # add room id to contact object to make multigraph more 
    # explicit  

    rooms_d = rooms.Day_to_Rooms[d]

    #id_to_contacts = N_lists_in.id_to_contacts

    #iterate through rooms 
     # TODO: update input arguments to take a list of rooms to update. 
    for room_id in room_ids_to_update#(room_id, room) in Rooms_d
        # iterate through agents in room
        room = rooms_d[room_id] 
        for source_id in room.agent_ids
            # assign out edges to all other agents 
            #(assign each a weight of 1 b/c it's a multigraph)



            # note: when adapting this to update existing edge lists, we have to first clear all old contacts FROM this room_id
            # use conditional indexing 
            contacts_to_keep = (c->c.room_id).(N_lists_in.id_to_contacts[source_id]) .!= room_id
            N_lists_in.id_to_contacts[source_id] = N_lists_in.id_to_contacts[source_id][contacts_to_keep] 
            # this works for the agents who have not been removed, but does not update 
            # the network for the removed agent (because they're not in the room to iterate over)

            # this will modify the contacts to clear all those with the current room_id
            # allowing us to then re-impute them from the current room assignments. 

            for target_id in room.agent_ids
                if target_id != source_id
                    weight = 1.0
                    contact = Networks_RACF.Contact(target_id, d, weight)
                    contact.room_id = room_id 
                    if haskey(N_lists_in.id_to_contacts, source_id)
                        push!(N_lists_in.id_to_contacts[source_id], contact) 
                    else
                        N_lists_in.id_to_contacts[source_id] = Array{Networks_RACF.Contact_T, 1}()
                        push!(N_lists_in.id_to_contacts[source_id], contact)
                    end 

                end
            end
        end
    end
end


end