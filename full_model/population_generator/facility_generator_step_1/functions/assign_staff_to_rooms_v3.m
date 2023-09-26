function [workers_out, rooms_out] = ...
                        assign_staff_to_rooms_v3(worker_struct,...
                                                 rooms,...
                                                 staff_per_room )

% CZ: 2022 09 01, sometimes, when there are many staff, it's possible to
% have staff who don't get assigned to rooms. fixed by checking at the end.
% 

n_rooms = size(rooms, 2);

room_ids = 1:n_rooms;


for d = 1:7
    
    workers_tbl = struct2table(worker_struct);
    
    all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);
    
    %zero room numbers for today, for fair assignment
    all_workers_today.n_rooms = zeros(size(all_workers_today.n_rooms));
    
    room_ids_rand = room_ids(randperm(size(rooms,2)));
    %scramble room ids for iterative assignmnet
    
    rm_id_list = room_ids_rand;

    for rm_i = rm_id_list %randomly iterate through rooms
         
        ids_working_r_today = [];
        
        for i = 1:staff_per_room
            
            %determine set of selectable ids for new staff assignment
            selectable_ids = [];
            
            % find staff working fewest rooms
            min_workload = min(all_workers_today.n_rooms);
            
            while isempty(selectable_ids)
                
                % preference anyone tied for lowest workload
                ids_min_workload = ...
                   all_workers_today.id(all_workers_today.n_rooms ==...
                                                            min_workload);
                
                % list of staff who have worked the room on previous days
                check_who_worked_before = @(x) is_room_there(x, rm_i, d);
                
                inds = cellfun(check_who_worked_before,...
                        all_workers_today.rooms, 'UniformOutput', false);
                    
                ids_worked_r_before = all_workers_today.id( [inds{:}] == 1);
                
                % preference those who worked the room before
                % but are not working the room today
                if isempty(ids_worked_r_before)
                    % if nobody worked the room before, then pick from
                    % anyone who is not working the room today
                    all_ids_today = all_workers_today.id;
                    available_ids = all_ids_today(~ismember(all_ids_today,...
                                                  ids_working_r_today));
                else
                    available_ids =...
                        ids_worked_r_before(~ismember(ids_worked_r_before,...
                                                      ids_working_r_today));
                end
                
                %enforces workload constraint
                selectable_ids = available_ids(ismember(available_ids,...
                                               ids_min_workload));
                
                %relaxes constraint that staff member worked the room
                %before (still selectable as long as their workload is low)
                if isempty(selectable_ids)
                    selectable_ids = ...
                        ids_min_workload(~ismember(ids_min_workload,...
                                                   ids_working_r_today));
                end
                
                % if all staff with minimal workload are already assigned,
                % increase the acceptable workload in order to assign 
                if isempty(selectable_ids)
                    min_workload = min_workload + 1;
                end
            end
            
            selected_worker = ...
                selectable_ids(randperm(numel(selectable_ids), 1)); 
            % random selection
            
            ids_working_r_today = [ids_working_r_today, selected_worker];
            
            %%%%%
            
            %assign selected staff member to room:
         
            rooms(rm_i).workers_general = [rooms(rm_i).workers_general;...
                                           [d, selected_worker, true]];
         
            worker_struct(selected_worker).rooms =...
                [worker_struct(selected_worker).rooms;...
                [d, rooms(rm_i).id, true] ];
            
            worker_struct(selected_worker).n_rooms =...
                worker_struct(selected_worker).n_rooms + 1 ;
            
            all_workers_today.n_rooms(all_workers_today.id == selected_worker) = ...
                all_workers_today.n_rooms(all_workers_today.id == selected_worker) + 1;
            
        end
    end
    
    % 2022 09 01
    % ensure that every agent present is assigned at least one room: 
    % note that this means rooms MAY have more than the minimum number
    % assigned. 
    
    %update worker structure: 
    workers_tbl = struct2table(worker_struct);
    all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);
    
    n_workers_today = size(all_workers_today, 1);
    
    % check if worker has been assigned a room:
    for i = 1:n_workers_today
        selected_worker = all_workers_today.id(i);
        rooms_i = all_workers_today.rooms{i};
        
        no_rooms_assigned_flag = false;
        if isempty(rooms_i)
            no_rooms_assigned_flag = true;
        elseif isempty(rooms_i(rooms_i(:, 1) == d, :))
            no_rooms_assigned_flag = true;
        end
        
        
        if no_rooms_assigned_flag
            
            rm_i = randperm(n_rooms, 1);
            
            %assign selected general staff member to random room:
            rooms(rm_i).workers_general = [rooms(rm_i).workers_general; [d, selected_worker, true]];

            
           worker_struct(selected_worker).rooms =...
                [worker_struct(selected_worker).rooms; [d, rooms(rm_i).id, true] ];
            
            worker_struct(selected_worker).n_rooms =...
                worker_struct(selected_worker).n_rooms + 1 ;
            
        end
    
    end
       
        
        
   
    
    
end



workers_out = worker_struct;
rooms_out = rooms;

end

%% helper functions

function a = is_room_there(x, id, d)
if ~isempty(x)
    rooms = x(:, 2);
    days = x(:, 1);
    rooms_to_check = rooms(days < d);
    a = ismember(id, rooms_to_check);
else
    a = 0;
end
end

