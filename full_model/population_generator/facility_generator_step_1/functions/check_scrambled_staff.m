% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function [H_out] = check_scrambled_staff(staff_ordered, staff_scrambled)


n_different_oo = 0;
n_tot_oo = 0;

n_different_os = 0;
n_tot_os = 0;

n = numel(staff_ordered);

for i = 1:n
    
    rooms_vs_d_ordered = {};
    rooms_vs_d_scrambled = {};
    
    iter = 0;
    
    for d = 1:7
        if staff_ordered(i).roster(d)
            iter = iter + 1;
            rooms_d_ordered = sort(staff_ordered(i).rooms(staff_ordered(i).rooms(:, 1) == d, 2));
            rooms_d_scrambled = sort(staff_scrambled(i).rooms(staff_ordered(i).rooms(:, 1) == d, 2));
            rooms_vs_d_ordered{iter} = rooms_d_ordered;
            rooms_vs_d_scrambled{iter} = rooms_d_scrambled;
            
            
        end
    end
    
    [n_tot_oo_i, n_unique_oo_i] = compare_room_sets(rooms_vs_d_ordered, rooms_vs_d_ordered, iter) ;
    
     n_different_oo = n_different_oo + n_unique_oo_i;
     n_tot_oo = n_tot_oo + n_tot_oo_i;
     
     [n_tot_os_i, n_unique_os_i] = compare_room_sets(rooms_vs_d_ordered, rooms_vs_d_scrambled, iter) ;
    
     n_different_os = n_different_os + n_unique_os_i;
     n_tot_os = n_tot_os + n_tot_os_i;
    
    
end

% H_out = (n_different_os - n_different_oo) / n_tot_os;
H_out = n_different_os / n_tot_os;


end

function [n_tot, n_different] = compare_room_sets(rooms_1, rooms_2, max_iter)

n_different = 0;
n_tot = 0;

for d_1 = 1:max_iter
    
    r_1 = rooms_1{d_1};
    
    for d_2 = 1:max_iter
        
        r_2 = rooms_2{d_2};
        
        if d_1 ~= d_2
        n_tot = n_tot + numel([r_1; r_2]);
        
        n_different = n_different + sum(~ismember(r_1, r_2)) +  sum(~ismember(r_2, r_1));
        end
        
    end
end

end
