% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function [residents_out, staff_med_out, staff_gen_out, rooms_out] = generate_agents_hyp_fac_v0(n_res, n_staff, facility_info, seed)

rng(seed)

%adapted from facility generator, v07

n = n_res;
k = n_staff;

P_med = 0.30; %aged care census nurses + allied care
k_med = ceil(P_med * k);
k_gen = k - k_med;

%these are guesses...
k_min_gen = ceil(k_gen / 5);
k_min_med = ceil(k_med / 5);

P_shared = facility_info.prop_res_shared_rooms;

P_1 = 1 - P_shared;
P_2 = P_shared;%0.3;
P_3 = 0;
P_4 = 0;%0.2;
%must add to 1, and P_1 always = P_shared - 1.
% must also be self-consistent with resident population
% NOTE: don't have much info, so we'll just have 1 or 2-person rooms for
% simplicity

P_h = facility_info.prop_residents_high_chc_needs;
P_b = facility_info.prop_residents_high_beh_needs;
P_a = facility_info.prop_residents_high_adl_needs;

% probability that residents in shared rooms are assigned to rooms
% with other residents with same needs.
%NOTE: little info, keeping at 0.5 for now.
P_same = 0.5;

%check aged 2020 care census:
% table 2.2
% p_FT = [50 / [50 + 144]
% p_PT = 144 / [50 + 144]
% assert: p_FT + p_3 + p_2 = 1
% gives p_2 = 0.38
% p_3 = 0.36
roster = [[1, 2, 3, 4, 5, 6, 7];...
    [0, 0.38, 0.36, 0, 0.26, 0, 0]];


% how many staff service each room.
gen_staff_per_room = 4;
med_staff_per_room = 1;

%scrambling levels (leave at 0, roster is messy enough)
H_med = 0.0 ;
H_gen = 0.0 ;

%% produce a unique label for outputs
% date_tag = ['tstamp_' datestr(datetime, 'YYYY_mm_dd_HH_MM_SS')];
% seed_tag = ['seed_' num2str(seed)];
% label = [version_tag, '_' seed_tag '_' date_tag];
%
% if save_config_flag
%     save(['config_MAT_' label]);
% end




%% TODO: output a table containing configurational details


%% STEP 1 initialise residents and workers
residents = struct;
workers_general = struct;
workers_medical = struct;
for i = 1:n
    residents(i).id = i;
    residents(i).chc_needs = false;
    residents(i).beh_needs = false;
    residents(i).adl_needs = false;
    residents(i).room_id = NaN;
end

%% STEP 1a generate work schedule from roster (constrained by k_min)
[schedule_general, success_flag_gen] = schedule_from_roster(k_gen, k_min_gen, roster);
[schedule_medical, success_flag_med] = schedule_from_roster(k_med, k_min_med, roster);

if ~success_flag_gen || ~success_flag_med
    return
end


for i = 1:k_gen
    workers_general(i).id = i;
    workers_general(i).medical = false;
    workers_general(i).roster = schedule_general(i, :);
    workers_general(i).rooms = [];
    workers_general(i).n_rooms = 0;
end

for i = 1:k_med
    workers_medical(i).id = i;
    workers_medical(i).medical = true;
    workers_medical(i).roster = schedule_medical(i, :);
    workers_medical(i).rooms = [];
    workers_medical(i).n_rooms = 0;
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
%P_same = 0.5;

%STEP 3a: create rooms

P_m = [P_1, P_2, P_3, P_4];
%counts of room by size (index corresponds to number of residents occupying room)

% CZ: 2022 08 22: need to round proportions for integer capacity: 

m_i = zeros(size(P_m)) ;

cap = 0;

for i = 1:size(P_m, 2)
    
    % round up for odd and down for even
    if ~mod(i, 2)
        m_i(i) = floor(P_m(i) * n / i);
    else
        m_i(i) = ceil(P_m(i) * n / i);
    end
    cap = cap + (m_i(i) * i);
end

%CZ: 2022 09 01 - rounding fails to produce enough rooms sometimes
% if there are not enough rooms, generate one extra: 
if cap < n
    m_i(1) = m_i(1) + 1;
end


prop_same = Inf;

% CZ: 2022 09 28, disabling this constraint for the generation of
% hypothetical facilities. 
tolerance = 1.0; % setting tolerance to 1 means any configuration will pass. 

discrepancy = Inf;
while discrepancy >= tolerance
    
    %create rooms - TODO: this for loop can be moved outside the while loop, and replaced
    %with a reinitialisation of the resident lists for each room.
    rooms = struct;
    room_id = 0;
    
    % 2022 08 22 : check that capacity is equal to or greater than the
    % number of residents: 
    n_check = 0;
    
    for sz = 1:size(m_i, 2)
        n_sz = m_i(sz) ;
        for r = 1:n_sz
            room_id = room_id + 1;
            rooms(room_id).id = room_id;
            rooms(room_id).capacity = sz;
            rooms(room_id).residents = [];
            rooms(room_id).workers_general = [];
            rooms(room_id).workers_medical = [];
            rooms(room_id).chc_needs = 0;
            n_check = n_check + sz;
        end
    end
    
        
    % STEP 3b assign residents to rooms
    
    % iterate through rooms, populate to capacity selecting from
    % residents with similar needs with probability p_same
    
    residents_placed = false(n, 1);
    room_ids = 1:size(rooms,2);
    room_ids_rand = room_ids(randperm(size(rooms,2)));
    
    % 2022 08 22 modifying for robustness: 
    % keep track of number placed, and make sure all residents get rooms. 
    
    for rm_i = room_ids_rand % shuffle room indices so that larger rooms are not always filled last
        % select the first resident
        r_i = 1;
        while residents_placed(r_i)
            r_i = randperm(n, 1);
        end
        residents_placed(r_i) = true;
        residents(r_i).room_id = rooms(rm_i).id;
        rooms(rm_i).residents = [rooms(rm_i).residents, residents(r_i).id];
        rooms(rm_i).capacity = rooms(rm_i).capacity - 1;
        
        % now place additional residents subject to the probability of grouping
        % individuals with similar needs, P_same
        r_j = r_i;
        while rooms(rm_i).capacity > 0
            %evaluate P_same
            pool = residents(~residents_placed);
            if rand() < P_same
                %NOTE: for now, confining this to medical needs, but could
                %expand to other types of needs.
                pool_same = pool([pool.chc_needs] == residents(r_j).chc_needs);
                
                % if there are no similar residents, then pick anyone who is unassigned
                if ~isempty(pool_same)
                    pool = pool_same;
                end
                
            else
                pool_diff = pool([pool.chc_needs] ~= residents(r_j).chc_needs);
                
                % if there are no different residents, then pick anyone who is unassigned
                if ~isempty(pool_diff)
                    pool = pool_diff;
                end
                
            end
            
            while residents_placed(r_j)
                ids = [pool.id];
                r_j = ids(randperm(numel(pool), 1));
            end
            
            residents_placed(r_j) = true;
            residents(r_j).room_id = rooms(rm_i).id;
            rooms(rm_i).residents = [rooms(rm_i).residents, residents(r_j).id];
            rooms(rm_i).capacity = rooms(rm_i).capacity - 1;
            
        end
    end
    
    %% TEST 1: check that rooms were assigned according to P_same
    prop_same = test_room_assignment(P_same, rooms, residents);
    % not perefect, but should be better with larger numbers etc.
    % while loop enforces a pre-defined tolerance level.
    discrepancy = abs(P_same - prop_same) * P_shared;
    disp(['discrepancy: abs(P_same - prop_same) x P_shared = ' num2str(abs(prop_same - P_same) * P_shared)])
end

%Step 3b. Classify rooms based on needs levels
for rm_i = 1:size(rooms, 2)
    
    room_i = rooms(rm_i);
    residents_i = rooms(rm_i).residents;
    
    if sum([residents(residents_i).chc_needs])
        rooms(rm_i).chc_needs = rooms(rm_i).chc_needs + 1;
    end
    
end

%% STEP 4: assign workers to rooms (note rostering is already done in 1a)
% each room needs workers based on needs
% assignment of workers to rooms may vary on each day of the schedule

random_rm_order_flag = 0;
%false seems to give more consistency (need to check that it doesn't create
%network artefacts though)

% assign general staff.
[workers_general, rooms] = ...
    assign_staff_to_rooms(workers_general, rooms, gen_staff_per_room, 'general', random_rm_order_flag);

% assign medical staff
[workers_medical, rooms] = ...
    assign_staff_to_rooms(workers_medical, rooms, med_staff_per_room, 'medical', random_rm_order_flag);




%% STEP 5: scramble workers randomly based on H_staff parameter.
% the parameter H_staff in [0, 1] represents the level of consistency in staff
% assignment from day to day, 0 is maximal regularity and 1 is random
% H_med = 1.0;
% H_gen = 1.0;

if H_med > 0 || H_gen > 0
    
    [workers_general_scrambled, rooms_scrambled] = scramble_room_assignments_v3(workers_general, rooms, H_gen, 'general');
    [workers_medical_scrambled, rooms_scrambled] = scramble_room_assignments_v3(workers_medical, rooms_scrambled, H_med, 'medical');
    
    % test scrambled configuration for quality control.
    H_out_general = check_scrambled_staff(workers_general, workers_general_scrambled);
    H_out_medical = check_scrambled_staff(workers_medical, workers_medical_scrambled);
    
    residents_out = residents;
    staff_gen_out = workers_general_scrambled;
    staff_med_out = workers_medical_scrambled;
    rooms_out = rooms_scrambled;
    
else
    
    residents_out = residents;
    staff_gen_out = workers_general;
    staff_med_out = workers_medical;
    rooms_out = rooms;
    
end

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