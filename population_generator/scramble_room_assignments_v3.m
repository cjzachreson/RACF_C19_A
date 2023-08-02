function [staff_out, rooms_out] = scramble_room_assignments_v3(staff, rooms, H, staff_type)

%2022 08 05, empty field 'workers' removed from rooms output struct. 

staff_out = staff;
rooms_out = rooms;

general_staff_flag = strcmp(staff_type, 'general');
medical_staff_flag = strcmp(staff_type, 'medical');
if ~(general_staff_flag || medical_staff_flag)
    disp('warning: invalid staff type staff not assigned to rooms.')
    return
end

% idea is to randomise the regularity of staff <-> room assignments

%iterate through rooms and swap staff with random staff from other rooms

room_ids = 1:size(rooms,2);

assignments = [];

for d = 1:7
    
    % H defines the number of scramblings
    assignments_d = assignment_list(rooms, d, general_staff_flag, medical_staff_flag);
    
    n_assignments = size(assignments_d, 1);
    
    n_to_scramble_today = floor(H * n_assignments / 2) * 2;
    
    n_scrambled_today = 0;
    
    w_i = NaN;
    w_j = NaN;
    rm_i = NaN;
    rm_j = NaN;
    
    % randomise the assignment list
    assignments_d_rand = assignments_d(randperm(size(assignments_d, 1)), :);
    scrambled_assignments_d = [];
    
    
    max_iter = n_to_scramble_today * 100;
    
    n_iter = 0;
    
    
    while n_scrambled_today < n_to_scramble_today 
        
        n_iter = n_iter + 1;
        if n_iter == max_iter
            %n_scrambled_today
            % reset
            assignments_d_rand = assignments_d(randperm(size(assignments_d, 1)), :);
            scrambled_assignments_d = [];
            n_scrambled_today = 0;
            
            w_i = NaN;
            w_j = NaN;
            rm_i = NaN;
            rm_j = NaN;
            n_iter = 0;
            disp(['reached maximum iterations, trying again']);
        end
        
        
    
        i = randperm(size(assignments_d_rand, 1), 1);
        j = randperm(size(assignments_d_rand, 1), 1);
        
        
        w_i = assignments_d_rand(i, 3);
        rm_i = assignments_d_rand(i, 1);
        
        w_j = assignments_d_rand(j, 3);
        rm_j = assignments_d_rand(j, 1);
        
        accept_flag = check_exchange(d, w_i, w_j, rm_i, rm_j,...
                                     scrambled_assignments_d);%,...
                                     %assignments_d_rand);
        
        if accept_flag
            
            w_k = w_i;
            rm_k = rm_j;
            
            w_l = w_j;
            rm_l = rm_i;
            
            a_k = [rm_k, d, w_k, false];
            a_l = [rm_l, d, w_l, false];
            
            % add scrabmled assignment for i and j to the list
            scrambled_assignments_d = [scrambled_assignments_d; [a_k;a_l]]; 
            % remove old asignment for i and j from the list
            assignments_d_rand([i, j], :) = [];
            
            n_scrambled_today = n_scrambled_today + 2;
        
        end
        
    end
    
    %put all assignments back together and reset output stucts
    assignments_d = [scrambled_assignments_d; assignments_d_rand] ;
    
    assignments = [assignments; assignments_d];
    
    disp(['scrambling day ' num2str(d) ', scrambled ' num2str(n_scrambled_today) ' room assignments'])
    
end

[rooms_out, staff_out] = reset_assignments(general_staff_flag, medical_staff_flag, rooms_out, staff_out, assignments);



end

function w = worker_ids_r_d(worker_list_r, d)
workers_r_d = worker_list_r(worker_list_r(:, 1) == d, :);
not_scrambled = workers_r_d(:, 3);
w = workers_r_d(not_scrambled == 1, 2);
end

function f = check_exchange(d, w_i, w_j, rm_i, rm_j, scrambled_assignments)%, old_assignments)
f = false;
% for an exchange to be accepted, it must not put someone to a room they've
% already been assigned to
if w_i == w_j
    return
elseif rm_i == rm_j
    return
end
a_k = [rm_j, d, w_i, false];
a_l = [rm_i, d, w_j, false];

if ~isempty(scrambled_assignments)
    if ismember(a_k, scrambled_assignments, 'row')  || ismember(a_l, scrambled_assignments, 'row')
        return
    end
end

% it's OK for the new one to be in the old set, but not for the new one to
% be in the new set already. 

% if ismember(a_k, old_assignments, 'row')  || ismember(a_l, old_assignments, 'row')
%     return
% end

f = true;


end



function assignment_list = assignment_list(rooms, d, general_flag, medical_flag)

assignment_list = [];

for r = 1:numel(rooms)
    
    if general_flag
        worker_list_r = rooms(r).workers_general;
    elseif medical_flag
        worker_list_r = rooms(r).workers_medical;
    end
    
    workers_r_d = worker_list_r(worker_list_r(:, 1) == d, :);
    room_ids = repmat(r, [size(workers_r_d, 1), 1]);
    
    assignment_list = [assignment_list; [room_ids, workers_r_d]];
    
    
end
end


function [rooms_out, staff_out] = reset_assignments(general_staff_flag, ...
                                                    medical_staff_flag, ...
                                                    rooms_out, staff_out,...
                                                    assignments)


    % sort by day.
    assignments = sortrows(assignments, [2, 1, 3]) ;
                                                
    % clear all existing assignments
    % modified on 2022 08 05
    for r = 1:numel(rooms_out)
        if general_staff_flag
            rooms_out(r).workers_general = [];
        elseif medical_staff_flag
            rooms_out(r).workers_medical = [];
        end
    end
    
    for s = 1:numel(staff_out)
        staff_out(s).rooms = [];
    end

    for i = 1:size(assignments, 1)
        
        a_i = assignments(i, :);
        
        w_i = a_i(3);
        rm_i = a_i(1);
        d_i = a_i(2);
        
        staff_out(w_i).rooms = [staff_out(w_i).rooms; [d_i, rm_i, false]];
       

        if general_staff_flag
            
            rooms_out(rm_i).workers_general = [rooms_out(rm_i).workers_general; [d_i, w_i, false]];
            
        elseif medical_staff_flag
            
            rooms_out(rm_i).workers_medical = [rooms_out(rm_i).workers_medical; [d_i, w_i, false]];
            
        end
        
    end

end
