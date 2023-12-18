
clear all
close all

% specify input and output locations

input_date_label = 'release_test\2023_08_02'
output_date_label = 'release_test\2023_12_18'

current_dir = pwd();

%locate directory with synthetic population info: 
cd('../')
pop_gen_dir = pwd();
cd('./prepare_agent_files_step_2')

input_dirname = [pop_gen_dir,'\facility_generator_step_1\output\',...
                 'agents_Hypothetical_facilities\' input_date_label];

output_dirname = [pwd() '\output\agents_model_ready\' output_date_label];

fac_list_dirname = [pop_gen_dir, ...
                    '\facility_generator_step_1\input\facility_characteristics'];
                
fac_list_fname = [fac_list_dirname, '\hypothetical_facility_characteristics.csv'];
fac_list = readtable(fac_list_fname);

n_fac = size(fac_list, 1);


for i = 1:n_fac
    
    %linelist entry for this outbreak: 
    fac_i = fac_list(i, :);
   
    facility_label = [num2str(fac_i.service_id) '_hyp'];
    
    fac_path = ['\' facility_label];
    
    fac_input_dirname = [input_dirname, '\' fac_path];
    
    
    if ~isfolder(fac_input_dirname)
        disp(['NO AGENT INFO FOR OUTBREAK: ' num2str(i)])
        continue
    end
    
    
    fname_rooms = [fac_input_dirname '\', 'rooms.csv'];
    fname_workers_G = [fac_input_dirname '\', 'staff_g_facility.csv'];
    fname_workers_M = [fac_input_dirname '\', 'staff_m_facility.csv'];
    fname_residents = [fac_input_dirname '\', 'res_facility.csv'];
    
    outdir_fac_i = [output_dirname, '\facID_' facility_label];
    
    if ~isfolder(outdir_fac_i)
        mkdir(outdir_fac_i)
    end
    %%
    
    residents = residents_from_csv(fname_residents);
    workers_G = workers_from_csv(fname_workers_G);
    workers_M = workers_from_csv(fname_workers_M);
    
    rooms = rooms_from_csv(fname_rooms);
    
    n_rooms = size(rooms, 2);
    
    uniform_weight = 1;
    
    n_residents = size(residents, 2);
    n_workers_G = size(workers_G, 2);
    n_workers_M = size(workers_M, 2);
    
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
    
    
    workers_G_ids = [];
    for r = 1:n_rooms
        workers_g_r = rooms(r).workers_general(:, 2) + n_residents;
        %alter the ids of the general workers
        rooms(r).workers_general(:, 2) = workers_g_r;
        workers_G_ids = [workers_G_ids; workers_g_r];
    end
    workers_G_ids = sort(unique(workers_G_ids));
    workers_G_ids =  num2cell(workers_G_ids');
    [workers_G.id] = workers_G_ids {:};
    
    workers_M_ids = [];
    for r = 1:n_rooms
        workers_m_r = rooms(r).workers_medical(:, 2) + n_residents + n_workers_G;
        %alter the ids of the medical workers:
        rooms(r).workers_medical(:, 2) = workers_m_r;
        workers_M_ids = [workers_M_ids; workers_m_r];
    end
    
    workers_M_ids = sort(unique(workers_M_ids));
    workers_M_ids =  num2cell(workers_M_ids');
   
    [workers_M.id] = workers_M_ids {:};

    workers_G_full = workers_G;
    workers_M_full = workers_M;
    residents_full = residents;
    
    
    %%
    write_agents_to_csv_Julia(residents_full, [outdir_fac_i '\residents_for_Julia_test']);
    write_agents_to_csv_Julia(workers_G_full, [outdir_fac_i '\workers_G_for_Julia_test']);
    write_agents_to_csv_Julia(workers_M_full, [outdir_fac_i '\workers_M_for_Julia_test']);
    
    
    
end






