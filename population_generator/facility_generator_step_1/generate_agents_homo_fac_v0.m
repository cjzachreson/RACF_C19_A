function [residents_out, rooms_out] =...
    generate_agents_homo_fac_v0(n_res, facility_info, seed)

rng(seed)

%adapted from facility generator, v07

n = n_res;

P_h = facility_info.prop_residents_high_chc_needs;
P_b = facility_info.prop_residents_high_beh_needs;
P_a = facility_info.prop_residents_high_adl_needs;



%% TODO: output a table containing configurational details


%% STEP 1 initialise residents and workers
residents = struct;

for i = 1:n
    residents(i).id = i;
    residents(i).chc_needs = false;
    residents(i).beh_needs = false;
    residents(i).adl_needs = false;
    residents(i).room_id = NaN;
end


%% STEP 2: assign residents behavioural needs types
%deterministic assigment of proportion with high needs.
n_h = floor(P_h * n);
n_b = floor(P_b * n);
n_a = floor(P_a * n);

ids_h = randperm(n, n_h);
ids_b = randperm(n, n_b);
ids_a = randperm(n, n_a);

for i = 1:n_h
    id_i = ids_h(i);
    residents(id_i).chc_needs = true;
end

for i = 1:n_b
    id_i = ids_b(i);
    residents(id_i).beh_needs = true;
end

for i = 1:n_a
    id_i = ids_a(i);
    residents(id_i).adl_needs = true;
end

%% STEP 3: assign residents to rooms
% CZ: 2023 09 22 - assigning all residents to one big room 
% (homogeneous structure). 

rooms = struct;
room_id = 1;

rooms(room_id).id = room_id;
rooms(room_id).capacity = n;
rooms(room_id).residents = [];

% STEP 3b assign residents to rooms
for rm_i = 1
    % select the first resident
    for r_i = 1:n
        residents(r_i).room_id = rooms(rm_i).id;
        rooms(rm_i).residents = [rooms(rm_i).residents, residents(r_i).id];
        rooms(rm_i).capacity = rooms(rm_i).capacity - 1;
    end
end



%% STEP 4: assign workers to rooms (note rostering is already done in 1a)
% each room needs workers based on needs
% assignment of workers to rooms may vary on each day of the schedule
residents_out = residents;
rooms_out = rooms;
    
%test writing script
%
% workers_general_label = ['workers_general_' label];
% workers_medical_label = ['workers_medical_' label];
% residents_label = ['residents_' label];
% rooms_label = ['rooms_' label];

% write_rooms_to_csv(rooms_scrambled, rooms_label);
% write_workers_to_csv(workers_general_scrambled, workers_general_label);
% write_workers_to_csv(workers_medical_scrambled, workers_medical_label);
% write_residents_to_csv(residents, residents_label);

end