 

clear all
close all

% add function directory to path: 
f_dir = [pwd() '\functions'];
addpath(f_dir);

%set rng seed
seed = 1;

%% load in data

facility_data_fname = ['input\facility_characteristics\',...
                       'hypothetical_facility_characteristics_v3.csv'];
facility_info = readtable(facility_data_fname);

%%output locations:
date_label = '2023_09_27';

output_dir = [pwd(), '\output\agents_Hypothetical_facilities_v3'];


%% specify input sets:
n_facilities = size(facility_info, 1);

for i = 1:n_facilities
    
    %facility index
    fac_i = facility_info(i, :);
    
    fac_id = fac_i.id;

    fac_id_label = ['facID_' num2str(fac_id)];
     
    % read in facility info used in pop structure generator
    i_facility_info = facility_info(facility_info.id == fac_id, :);
    
    i_n_residents = fac_i.n_residents;
    i_n_staff = fac_i.n_staff;
    
    
    [i_residents_fac,...
     i_staff_gen_fac,...
     i_rooms_fac  ]  = ...
                        generate_agents_hyp_fac_v3(i_n_residents,...
                                                   i_n_staff,...
                                                   i_facility_info,...
                                                   seed);
    
    % add outbreak id and key to agent and room structures before writing (for final linking of data):
    n_res = size(i_residents_fac, 2);
    n_wrk_g = size(i_staff_gen_fac, 2);
    n_rooms = size(i_rooms_fac, 2);
    
    
    
    %% make a directory for the set of agents from the hypothetical facility:
    fac_i_agent_list_dirname = [num2str(fac_i.id), '_hyp'];
    out_dirname = [output_dir, '\' date_label '\' fac_i_agent_list_dirname];
    if ~isfolder(out_dirname)
        mkdir(out_dirname)
    end
    
    %% write resident, worker, and room datasets for each outbreak:
   
    res_fac_label = [out_dirname, '\res_facility'];
    wrk_fac_label = [out_dirname, '\staff_facility'];
    rm_fac_label = [out_dirname, '\rooms'];
    
    
    %note that these are multidimensional fields, so they require customised
    %read-in scripts (which are included in the network generator).
    
    %write workers to csv (general)
    write_workers_to_csv(i_staff_gen_fac, wrk_fac_label)
    
    %write residents to csv
    write_residents_to_csv(i_residents_fac, res_fac_label)
    
    %write rooms to csv
    write_rooms_to_csv(i_rooms_fac, rm_fac_label)
    
end

% format for julia read-in: 

% specify input and output locations
input_dirname = [output_dir, '\', date_label];

output_dirname = [output_dir, '\agents_model_ready_v3\' date_label];

format_agents_for_Julia_readin(input_dirname, output_dirname, facility_data_fname)


