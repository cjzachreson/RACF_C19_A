
function format_agents_for_Julia_readin(input_dirname, output_dirname, fac_list_fname)

fac_list = readtable(fac_list_fname);

n_fac = size(fac_list, 1);


for i = 1:n_fac
    
    %linelist entry for this outbreak: 
    fac_i = fac_list(i, :);
   
    facility_label = [num2str(fac_i.id) '_hyp'];
    
    fac_path = ['\' facility_label];
    
    fac_input_dirname = [input_dirname, '\' fac_path];
    
    
    if ~isfolder(fac_input_dirname)
        disp(['NO AGENT INFO FOR FACILITY: ' num2str(i)])
        continue
    end
    
    
    fname_rooms = [fac_input_dirname '\', 'rooms.csv'];
    fname_workers = [fac_input_dirname '\', 'staff_facility.csv'];
    fname_residents = [fac_input_dirname '\', 'res_facility.csv'];
   
    %fname_residents_imstat = [fac_input_dirname '\', 'res_imstat.csv'];
    %fname_workers_imstat = [fac_input_dirname '\', 'staff_imstat.csv'];
    
    
    outdir_fac_i = [output_dirname, '\facID_' facility_label];
    
    if ~isfolder(outdir_fac_i)
        mkdir(outdir_fac_i)
    end
    %%
    
    residents = residents_from_csv(fname_residents);
    workers = workers_from_csv(fname_workers);
    
    rooms = rooms_from_csv(fname_rooms);
    
    n_rooms = size(rooms, 2);
    
    n_residents = size(residents, 2);
    
    % compile lists of workers and residents from structures:
    
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
    
    
    workers_ids = [];
    for r = 1:n_rooms
        workers_r = rooms(r).workers_general(:, 2) + n_residents;
        %alter the ids of the general workers
        rooms(r).workers_general(:, 2) = workers_r;
        workers_ids = [workers_ids; workers_r];
    end
    workers_ids = sort(unique(workers_ids));
    workers_ids =  num2cell(workers_ids');
    [workers.id] = workers_ids {:};
    

    workers_full = workers;
    residents_full = residents;
    
    
    %%
    write_agents_to_csv_Julia(residents_full, [outdir_fac_i '\residents_for_Julia_test']);
    write_agents_to_csv_Julia(workers_full, [outdir_fac_i '\workers_for_Julia_test']);
       
    
end

end






