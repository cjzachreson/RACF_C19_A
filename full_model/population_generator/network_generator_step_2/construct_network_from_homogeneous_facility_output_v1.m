
clear all
close all

% specify input and output locations

input_date_label = '2023_09_22'
output_date_label = '2023_09_22'

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


%% Helper functions
function E_keys = STD_2_unique_keys(E_list, sep)

%separator
E_keys = {};
for e_i = 1:size(E_list, 1)
    
    e = E_list(e_i, :);
    
    key_i = [num2str(e(1)), sep, num2str(e(2)), sep, num2str(e(3))];
    
    E_keys{e_i, 1} = key_i;
    
end

E_keys = unique(E_keys);

end

function STD_W_map = Multigraph_2_weightmap(STD_keys, E_list_mgraph, sep)

STD_W_map = containers.Map(STD_keys, zeros(size(STD_keys, 1), 1));

for e_i = 1:size(E_list_mgraph, 1)
    
    e = E_list_mgraph(e_i, :);
    
    key_i = [num2str(e(1)), sep, num2str(e(2)), sep, num2str(e(3))];
    
    weight_e_i = e(4);
    
    STD_W_map(key_i) = STD_W_map(key_i) + weight_e_i;
    
end

end

function W_E_list = STD_weightmap_2_Elist(STD_Wmap, sep)

key_set = keys(STD_Wmap)';

W_E_list = NaN(size(key_set, 1), 4);

for k_i = 1:size(key_set, 1)
    
    key_i = key_set{k_i};
    
    w_i = STD_Wmap(key_i);
    
    STD = strrep(key_i, sep, ',');
    
    ST_num = eval([ '[' STD ']' ]);
    
    W_E_list(k_i, :) = [ST_num, w_i];
    
end

W_E_list = sortrows(W_E_list, [1,3,2]);

end

function NList = WNet_2_NList(WE_list)

S_list = unique(WE_list(:, 1));

NList = cell([numel(S_list), 2]);

for s_i = 1:numel(S_list)
    
    source_id = S_list(s_i);
    NList{s_i, 1} = source_id;
    
    target_ids = WE_list( WE_list(:, 1) == source_id, 2);
    target_days = WE_list(WE_list(:, 1) == source_id, 3);
    target_weights = WE_list(WE_list(:, 1) == source_id, 4);
    NList{s_i, 2} = [target_ids, target_days, target_weights];
    
end

end


function struct_out = struct_cat(s1, s2)

c1 = struct2cell(s1);
c2 = struct2cell(s2);

fields_s1 = fieldnames(s1);
fields_s2 = fieldnames(s2);

fields_12 = {fields_s1{:}, fields_s2{:}}';

n = size(s1, 2);

for i = 1:n
    
    c12 = {c1{:, :, i}, c2{:, i}}';
    
    if i == 1
    
    struct_out = cell2struct(c12, fields_12);
    
    else
        
        struct_out(i) = cell2struct(c12, fields_12);
        
    end

end
end




