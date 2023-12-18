
clear all
close all

% specify input and output locations

input_date_label = '2023_09_22'
output_date_label = '2023_12_16'

current_dir = pwd();

%locate directory with synthetic population info: 
cd('../')
pop_gen_dir = pwd();
cd('./network_generator_step_2')

input_dirname = [pop_gen_dir,'\facility_generator_step_1\',...
                 'agents_Homogeneous_test_pops\' input_date_label];

output_dirname = [pwd() '\agents_model_ready\' output_date_label];

fac_list_dirname = [pop_gen_dir, ...
                    '\facility_generator_step_1\facility_characteristics'];
fac_list_fname = [fac_list_dirname, '\homogeneous_facility_characteristics.csv'];
fac_list = readtable(fac_list_fname);

n_fac = size(fac_list, 1);


for i = 1:n_fac
    
    %linelist entry for this outbreak: 
    fac_i = fac_list(i, :);
   
    facility_label = [num2str(fac_i.service_id) '_homo'];
    
    fac_path = ['\' facility_label];
    
    fac_input_dirname = [input_dirname, '\' fac_path];
    
    
    if ~isfolder(fac_input_dirname)
        disp(['NO AGENT INFO FOR FACILITY: ' num2str(i)])
        continue
    end
    
    
    fname_rooms = [fac_input_dirname '\', 'rooms.csv'];
    fname_residents = [fac_input_dirname '\', 'res_facility.csv'];
   
    
    
    outdir_fac_i = [output_dirname, '\facID_' facility_label];
    
    if ~isfolder(outdir_fac_i)
        mkdir(outdir_fac_i)
    end
    %%
    
    residents = residents_from_csv(fname_residents);
    
    rooms = rooms_from_csv(fname_rooms);
    
    n_rooms = size(rooms, 2);
    
    uniform_weight = 1;
    
    %NOTE: for the network, medical workers and general workers need to come from
    % the same list of identifiers, this can be done by adding the total number of general
    %workers to the ids of the medical workers (or vice versa).
    % eg: ID_medical(new) = ID_medical(original) + N(workers_general);
    % need to do the same with residents... all IDs must be unique.
    
    n_residents = size(residents, 2);
    
    % compile lists of workers and residents from structures:
    % note, there is redundancey here, for clarity.
    
    resident_ids = [];
    for r = 1:n_rooms
        residents_r = rooms(r).residents;
        %no alteration of resident IDs... these are counted first.
        resident_ids = [resident_ids; residents_r'];
    end
    resident_ids = sort(unique(resident_ids));
    resident_ids =  num2cell(resident_ids');
    
    
    [residents.id] = resident_ids{:};
    %note that (for now) the above step is not necessary.
    % b.c. residents have only one room assignment each. Leaving in for
    % robustness.
    
   
    residents_full = residents;
    
    
    %%
    write_agents_to_csv_Julia(residents_full, [outdir_fac_i '\residents_for_Julia_test']);
    

    % disabling all network generators - these are now constructed in
    % Julia from the agent-room assignments. 

    
    
end



