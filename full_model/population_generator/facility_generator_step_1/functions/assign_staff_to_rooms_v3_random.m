function [workers_out, rooms_out] = ...
                        assign_staff_to_rooms_v3_random(worker_struct,...
                                                        rooms,...
                                                        min_staff_per_room,...
                                                        min_rooms_per_staff)



n_rooms = size(rooms, 2);

room_ids = 1:n_rooms;

% TODO: fix repeat assignment algorithm so it overwrites 
% appropriately. 

    for d = 1:7


        workers_tbl = struct2table(worker_struct);

        all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);

        n_workers_today = size(all_workers_today, 1);

        room_ids_rand = room_ids(randperm(size(rooms,2)));
        %scramble room ids for iterative assignmnet

        rm_id_list = room_ids_rand;

        for rm_i = rm_id_list %randomly iterate through rooms

            % select staff at random:
            idx = randperm(n_workers_today, min_staff_per_room);
            staff_rm_d = all_workers_today.id(idx);

            %assign selected staff member to room:

            for j = 1:min_staff_per_room

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

        % for agents assigned fewer than the minimum number of rooms, 
        % give additional assignments at random: 

        %update list of today's workers with new assignments: 
        workers_tbl = struct2table(worker_struct);
        all_workers_today = workers_tbl(workers_tbl.roster(:, d) == 1, :);
        for i = 1:n_workers_today

            staff_id_i = all_workers_today.id(i);

            rooms_i = all_workers_today.rooms{i};
            
            if isempty(rooms_i)
                rooms_i_d = [];
            else
                rooms_i_d = rooms_i(rooms_i(:, 1) == d);
            end
            
            n_rooms_i_d = size(rooms_i_d, 1);

            if n_rooms_i_d < min_rooms_per_staff
                n_add_rm = min_rooms_per_staff - n_rooms_i_d;
                idx = randperm(n_rooms, n_add_rm);

                for j = 1:n_add_rm
                    rm_id_j = idx(j);
                    rooms(rm_id_j).workers_general = [rooms(rm_id_j).workers_general;...
                                                     [d, staff_id_i, true]];

                    worker_struct(staff_id_i).rooms =...
                        [worker_struct(staff_id_i).rooms;...
                        [d, rooms(rm_id_j).id, true] ];

                    worker_struct(staff_id_i).n_rooms =...
                        worker_struct(staff_id_i).n_rooms + 1 ;

                    all_workers_today.n_rooms(all_workers_today.id == staff_id_i) = ...
                        all_workers_today.n_rooms(all_workers_today.id == staff_id_i) + 1;

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

