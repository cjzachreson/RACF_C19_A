 

clear all
close all

%set rng seed
seed = 1;

%% load in data
date_label = '2023_09_22';


facility_data_fname = ['facility_characteristics\homogeneous_facility_characteristics.csv'];

facility_info = readtable(facility_data_fname);

%%output locations:

output_dir = [pwd(), '\agents_Homogeneous_test_pops'];


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
   
    % read in facility info used in pop structure generator
    i_facility_info = facility_info(facility_info.service_id == i_service_id, :);
    
    i_num_residents = fac_i.n_residents;
    i_num_staff = fac_i.n_staff;
    
    
    [i_residents_fac, i_rooms_fac  ]  = ...
        generate_agents_homo_fac_v0(i_num_residents, i_facility_info, seed);
    
    % add outbreak id and key to agent and room structures before writing (for final linking of data):
    n_res = size(i_residents_fac, 2);
    n_rooms = size(i_rooms_fac, 2);
    
    for res_i = 1:n_res
        i_residents_fac(res_i).outbreak_key = 'xxxx';
    end
    
    for rm_i = 1:n_rooms
        i_rooms_fac(rm_i).outbreak_key = 'xxxx';
    end
    
    
    %% make a directory for the set of agents from the hypothetical facility:
    fac_i_agent_list_dirname = [num2str(fac_i.service_id), '_homo'];
    out_dirname = [output_dir, '\' date_label '\' fac_i_agent_list_dirname];
    if ~isfolder(out_dirname)
        mkdir(out_dirname)
    end
    
    %% write resident, worker, and room datasets for each outbreak:
   
    res_fac_label = [out_dirname, '\res_facility'];
    rm_fac_label = [out_dirname, '\rooms'];
    
    
    %note that these are multidimensional fields, so they require customised
    %read-in scripts (which are included in the network generator).
    
    %write residents to csv
    write_residents_to_csv(i_residents_fac, res_fac_label)
    
    %write rooms to csv
    write_rooms_to_csv(i_rooms_fac, rm_fac_label)
    
end

