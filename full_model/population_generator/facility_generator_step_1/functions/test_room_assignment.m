function [prop_same] = test_room_assignment(P_same, rooms, residents) 

% counts number of shared rooms assigned to residents with same needs

n_same = 0;
n_diff = 0;


for rm_i = 1:size(rooms,2)
    
    residents_rm = residents(rooms(rm_i).residents);
    if numel(residents_rm) > 1
        if sum(diff([residents_rm.chc_needs])) == 0
            n_same = n_same + 1;
        else
            n_diff = n_diff + 1;
        end
    end
    
end

n_tot = n_same + n_diff;
prop_same = n_same/n_tot;


        
    
    