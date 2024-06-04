% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function [schedule, success_flag] = schedule_from_roster(k, k_min, roster) 
success_flag = true;
max_attempts = 10000;

schedule = zeros(k, 7);

reject_flag = true;
n_attempts = 0;

while reject_flag 
    
    for i = 1:k
        n_days = randsample(roster(1, :), 1, true, roster(2,:));
        shift = randperm(7, 1);
        s_i = [ones(1, n_days), zeros(1, 7 - n_days)];
        s_i = circshift(s_i, shift);
        schedule(i, :) = s_i;
    end
   
    test = sum(schedule);
    if min(test) >= k_min
        reject_flag = false;
        disp('accepted a schedule') 
    end
    
    n_attempts = n_attempts + 1;
    
    if n_attempts == max_attempts
        disp('max attempts reached, trouble assigning schedule')
        reject_flag = false;
        success_flag = false;
    end
    
end
        
    
    