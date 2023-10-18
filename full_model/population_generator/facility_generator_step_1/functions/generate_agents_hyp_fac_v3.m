function [residents_out,...
          staff_out,...
          rooms_out] = ...
                        generate_agents_hyp_fac_v3(n_res,...
                                                   n_staff,...
                                                   facility_info,...
                                                   seed)

rng(seed)

n = n_res;
k = n_staff;

% minimum number of staff
% required on a given day 
k_min = ceil(k / 5);

% define work roster. 
% TODO: add proper reference. 
%check 2020 aged care census:
% table 2.2
% p_FT = [50 / [50 + 144]
% p_PT = 144 / [50 + 144]
% assert: p_FT + p_3 + p_2 = 1
% gives p_2 = 0.38
% p_3 = 0.36
roster = [[1, 2, 3, 4, 5, 6, 7];...
    [0, 0.38, 0.36, 0, 0.26, 0, 0]];


% best-guess. 
% how many staff service each room.
staff_per_room = facility_info.min_staff_per_room;
rooms_per_staff = facility_info.min_rooms_per_staff;
% TODO: check if this is per-day or per roster period 
% (I think it's per day). 

%scrambling levels
% for adding additional randomisation to the roster. 
% TODO: roster permutations should be made in the 
% simulator, not in the population generator. 
H = 0.0;

%% ***TODO: output a table containing configurational details


%% STEP 1 initialise residents and workers
residents = struct;
staff = struct;
for i = 1:n
    residents(i).id = i;
    residents(i).room_id = NaN;
end

%% STEP 1 generate work schedule from roster (constrained by k_min)
[staff_schedule,...
 success_flag_sched] = schedule_from_roster(k,...
                                            k_min,...
                                            roster);


if ~success_flag_sched
    return
end


for i = 1:k
    staff(i).id = i;
    staff(i).roster = staff_schedule(i, :);
    staff(i).rooms = [];
    staff(i).n_rooms = 0;
end


%% STEP 2: assign residents to rooms

%STEP 2a: create rooms

rooms = struct;

sz = facility_info.residents_per_room;

n_rooms = ceil(n / sz);

for r = 1:n_rooms
    room_id = r;
    rooms(room_id).id = room_id;
    rooms(room_id).capacity = sz;
    rooms(room_id).residents = [];
    rooms(room_id).workers_general = [];
end


% STEP 2b fill rooms with residents,
% subject to room capacity. 
res_id = 0;
rm_id = 1;
while res_id < n  

    % place residents 

    while rooms(rm_id).capacity > 0
        
        res_id = res_id + 1;
        
        if res_id > n
            break
        end
       
        residents(res_id).room_id = rooms(rm_id).id;
        rooms(rm_id).residents = [rooms(rm_id).residents, residents(res_id).id];
        rooms(rm_id).capacity = rooms(rm_id).capacity - 1;
    end
    
    rm_id = rm_id + 1;
    
end
    

%% STEP 4: assign workers to rooms (note rostering is already done in 1a)
% assignment of workers to rooms may vary on each day of the schedule

% assign general staff.
% note - this assignment system attempts to 
% distribute the workload among staff
% while keeping the room assignments consistent 
% from day to day (noting that this is not possible
% to do perfectly), optimal rostering is out of scope. 

% [staff, rooms] = ...
%                     assign_staff_to_rooms_v3(staff,...
%                                              rooms,...
%                                              staff_per_room);

[staff, rooms] = ...
                    assign_staff_to_rooms_v3_random(staff,...
                                                    rooms,...
                                                    staff_per_room,...
                                                    rooms_per_staff);
                                         
                                         





%% STEP 5: scramble workers randomly based on H_staff parameter.
% the parameter H_staff in [0, 1] represents the level of consistency in staff
% assignment to sets of rooms from day to day, 0 is maximal regularity and 1 is random
% H_med = 1.0;
% H_gen = 1.0;

if H > 0 
    
    [staff_scrambled, rooms_scrambled] = ...
            scramble_room_assignments_v3(staff,...
                                         rooms,...
                                         H);

    
    % test scrambled configuration for quality control.
    H_out = check_scrambled_staff(staff, staff_scrambled);
    
    residents_out = residents;
    staff_out = staff_scrambled;
    rooms_out = rooms_scrambled;
    
else
    
    residents_out = residents;
    staff_out = staff;
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