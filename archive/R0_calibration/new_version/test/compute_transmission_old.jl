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