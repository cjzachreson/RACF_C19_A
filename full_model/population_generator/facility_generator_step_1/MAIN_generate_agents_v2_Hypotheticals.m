 % Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

clear all
close all

addpath([pwd(), '\functions'])

%set rng seed
seed = 1;

%% load in data
date_label = '2023_08_02';


facility_data_fname = ['input\facility_characteristics\'...
                       'hypothetical_facility_characteristics.csv'];

facility_info = readtable(facility_data_fname);

%%output locations:

output_dir = [pwd(), '\output\agents_Hypothetical_facilities\release_test\'];


%% specify input sets:

%outbreak index to examine (can iterate through if desired):

% iterate through the first few outbreaks:

n_fac = size(facility_info, 1);

for i = 1:n_fac
    
    % line list entry for grouped outbreak (group corresponds to
    % classification of exemplar facility)
    fac_i = facility_info(i, :);
    
    i_service_id = fac_i.service_id;

    fac_id_label = ['facID_' num2str(i_service_id)];
    
    %placeholders for immunity status: 
    i_residents_immunity_status = table();
    
    i_staff_immunity_status = table();
   
    % read in facility info used in pop structure generator
    i_facility_info = facility_info(facility_info.service_id == i_service_id, :);
    
    i_num_residents = fac_i.n_residents;
    i_num_staff = fac_i.n_staff;
    
    
    [i_residents_fac, i_staff_med_fac, i_staff_gen_fac, i_rooms_fac  ]  = ...
        generate_agents_hyp_fac_v0(i_num_residents, i_num_staff, i_facility_info, seed);
    
    % add outbreak id and key to agent and room structures before writing (for final linking of data):
    n_res = size(i_residents_fac, 2);
    n_wrk_g = size(i_staff_gen_fac, 2);
    n_wrk_med = size(i_staff_med_fac, 2);
    n_rooms = size(i_rooms_fac, 2);
    
    for res_i = 1:n_res
        i_residents_fac(res_i).outbreak_key = 'xxxx';
    end
    
    for w_g = 1:n_wrk_g
        i_staff_gen_fac(w_g).outbreak_key = 'xxxx';
    end
    
    for w_m = 1:n_wrk_med
        i_staff_med_fac(w_m).outbreak_key = 'xxxx';
    end
    
    for rm_i = 1:n_rooms
        i_rooms_fac(rm_i).outbreak_key = 'xxxx';
    end
    
    
    %% make a directory for the set of agents from the hypothetical facility:
    fac_i_agent_list_dirname = [num2str(fac_i.service_id), '_hyp'];
    out_dirname = [output_dir, '\' date_label '\' fac_i_agent_list_dirname];
    if ~isfolder(out_dirname)
        mkdir(out_dirname)
    end
    
    %% write resident, worker, and room datasets for each outbreak:
   
    res_fac_label = [out_dirname, '\res_facility'];
    wrk_g_fac_label = [out_dirname, '\staff_g_facility'];
    wrk_m_fac_label = [out_dirname, '\staff_m_facility'];
    rm_fac_label = [out_dirname, '\rooms'];
    
    
    %note that these are multidimensional fields, so they require customised
    %read-in scripts (which are included in the network generator).
    
    %write workers to csv (general)
    write_workers_to_csv(i_staff_gen_fac, wrk_g_fac_label)
    
    %write workers to csv (medical)
    write_workers_to_csv(i_staff_med_fac, wrk_m_fac_label)
    
    %write residents to csv
    write_residents_to_csv(i_residents_fac, res_fac_label)
    
    %write rooms to csv
    write_rooms_to_csv(i_rooms_fac, rm_fac_label)
    
end

