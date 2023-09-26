function [workers_out, rooms_out] = ...
                        assign_staff_to_rooms_v3_random(worker_struct,...
                                                        rooms,...
                                                        staff_per_room )



n_rooms = size(rooms, 2);

room_ids = 1:n_rooms;

rooms_initial = rooms;
staff_initial = worker_struct;

% TODO: fix repeat assignment algorithm so it overwrites 
% appropriately. 

for d = 1:7
    
    repeat_assignment_flag = true
    
    while repeat_assignment_flag
    
        workers_tbl = struct2table(worker_struct);

        all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);

        n_workers_today = size(all_workers_today, 1);

        room_ids_rand = room_ids(randperm(size(rooms,2)));
        %scramble room ids for iterative assignmnet

        rm_id_list = room_ids_rand;

        for rm_i = rm_id_list %randomly iterate through rooms
            
            % select staff at random:
            idx = randperm(n_workers_today, staff_per_room);
            staff_rm_d = all_workers_today.id(idx);

            %assign selected staff member to room:

            for j = 1:staff_per_room

                staff_j = staff_rm_d(j);

                rooms(rm_i).workers_general = [rooms(rm_i).workers_general;...
                                               [d, staff_j, true]];

                worker_struct(staff_j).rooms =...
                    [worker_struct(staff_j).rooms;...
                    [d, rooms(rm_i).id, true] ];

                worker_struct(staff_j).n_rooms =...
                    worker_struct(staff_j).n_rooms + 1 ;

                all_workers_today.n_rooms(all_workers_today.id == staff_j) = ...
                    all_workers_today.n_rooms(all_workers_today.id == staff_j) + 1;

            end

        end



        % 2022 09 01
        % ensure that every agent present is assigned at least one room: 
        % note that this means rooms MAY have more than the minimum number
        % assigned. 


        repeat_assignment_flag = false; %reset to true below
        % if any workers were not assigned a room. 

        %update worker structure: 
        workers_tbl = struct2table(worker_struct);
        all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);


        % check if worker has been assigned a room:
        for i = 1:n_workers_today
           
            rooms_i = all_workers_today.rooms{i};

            no_rooms_assigned_flag = false;
            if isempty(rooms_i)
                no_rooms_assigned_flag = true;
            elseif isempty(rooms_i(rooms_i(:, 1) == d, :))
                no_rooms_assigned_flag = true;
            end


            if no_rooms_assigned_flag

                % repeat the sampling process
                repeat_assignment_flag = true;

            end

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

