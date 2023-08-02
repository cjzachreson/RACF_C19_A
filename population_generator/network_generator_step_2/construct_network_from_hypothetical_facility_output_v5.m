
clear all
close all

% specify input and output locations

input_date_label = '2022_09_28'
output_date_label = '2022_09_29'


pop_gen_dir = 'C:\Users\czachreson\Desktop\policy_work\RACF_2022\RACF_code\population_generator\';
%OB_agent_list_dir = [pop_gen_dir, input_date_label '\vaccination_timelines_step_2\outputs\'];

%just using this to get the list of outbreaks
%OB_resident_list_fname = [OB_agent_list_dir, 'OB_resident_immunity_status__rng_1_Tan_0.2_Tex_0.2_DtPro_2_DtMin_8_ver_v0.csv'];


%OB_residents = readtable(OB_resident_list_fname);

%output_dirname = [pwd() '\agents_OB_res_bkg_contact\' output_date_label]


input_dirname = [pop_gen_dir, '\' input_date_label,...
    '\facility_generator_step_3\agents_Hypothetical_facilities\' input_date_label];

output_dirname = [pwd() '\agents_model_ready_hypothetical\' output_date_label];

%iterate through outbreak set and construct networks

fac_list_dirname = [pop_gen_dir '\' input_date_label '\facility_classifier_step_0'];
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
   
    %fname_residents_imstat = [fac_input_dirname '\', 'res_imstat.csv'];
    %fname_workers_imstat = [fac_input_dirname '\', 'staff_imstat.csv'];
    
    
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
    
    %NOTE: for the network, medical workers and general workers need to come from
    % the same list of identifiers, this can be done by adding the total number of general
    %workers to the ids of the medical workers (or vice versa).
    % eg: ID_medical(new) = ID_medical(original) + N(workers_general);
    % need to do the same with residents... all IDs must be unique.
    
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
    
%     if size([residents.id], 2) ~= size(resident_ids, 2)
%         'check'
%     end
    
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
    
    if size([workers_M.id], 2) ~= size(workers_M_ids, 2)
        'check'
    end
    
    [workers_M.id] = workers_M_ids {:};
    
    %make a multigraph from the individuals who share rooms:
    
    %iterate through rooms and days and add an edge to the multigraph for each
    %individual sharing the space.
    
    %make a different network for each day:
    
    %edge_list = [];
    
    %CZ: in this version, the lists are concatenated with the day as an element
    %of each entry. (single network as output with time indicators on each edge).
    % this makes it easier to read in to julia for agent set construction.
    
    %Edge: [source, target, day, weight]
    
    % 2022 08 30 : add a 'super-room' containing all residents 
    % 2022 09 01 : removing supersource b/c this is implemented in the
    % simulation stage (saves a lot of memory) 
    
    %super_room_id = n_rooms + 1;
    %rooms(super_room_id).id = super_room_id;
    %rooms(super_room_id).capacity = 0;
    %rooms(super_room_id).residents = [resident_ids{:}]';
    %%for now, do not include workers in the super room. 
    %n_rooms = n_rooms + 1;
%     
%     for d = 1:7 %iterate through the days of the week:
%         
%         %iterate through rooms and create the multigraph for each day:
%         
%         for r = 1:n_rooms
%             
%             residents_r = rooms(r).residents';
%             
%             if ~isempty(rooms(r).workers_general)
%                 workers_gen_r_d = rooms(r).workers_general(rooms(r).workers_general(:, 1) == d, 2);
%             else
%                 workers_gen_r_d = [];
%             end
%             
%             if ~isempty(rooms(r).workers_medical)
%                 workers_med_r_d = rooms(r).workers_medical(rooms(r).workers_medical(:, 1) == d, 2);
%             else
%                 workers_med_r_d = [];
%             end
%             
%             
%             nodes_d_r = [residents_r; workers_gen_r_d; workers_med_r_d];
%             
%             edges_d_r = nchoosek(nodes_d_r, 2);
%             
%             % bi-directional edges, this is necessary because the neighbour
%             % lists will be used for filtering by source node.
%             edges_d_r = [edges_d_r ; fliplr(edges_d_r)] ;
%             
%             edges_d_r = [edges_d_r, ones(size(edges_d_r, 1), 1) * d, ones(size(edges_d_r, 1), 1) * uniform_weight];
%             
%             edge_list = [edge_list; edges_d_r];
%             
%         end
%         
%     end
    
    % convert edge list multigraphs to weighted sparse adjacency lists
    % (these will be exported as agent properties)
    
    %Multigraph_EL = sortrows(edge_list, [1, 3, 2]);
    %Weighted_Networks = {};
    %NeighbourLists = {};
    
    %first, make list of unique keys
    %sep = '#';
    
    %find set of unique edges in multigraph
    %STD_keys = STD_2_unique_keys(Multigraph_EL, sep);
    
    % accumulate multigraph edge list to weighted edge list:
    %STD_Wmap = Multigraph_2_weightmap(STD_keys, Multigraph_EL, sep);
    
    %convert back to S, T, W from weightmap
    %Weighted_Networks = STD_weightmap_2_Elist(STD_Wmap, sep);
    
    %neighbour lists will be formatted as: { [S] [T1; T2; T3...] }
    
    %convert to neighbour list
    % for speed in implementation of temporal networks, this will output
    % ALL of: multigraph edge lists; weighted edge lists; and neighbour lists.
    % TODO: determine if both mgraph and wgraph are needed (mgraph is
    % probably less general... given that a wgraph can be converted to
    % mgraph, the utility will depend on edge sampling routines).
    
    % adjacency list will be useful for conditional selection of edges
    % based on who is infected (i.e., excluding ST where S is not infected).
    
    %NeighbourLists = WNet_2_NList(Weighted_Networks);
    
    % print in a format easily readable into Julia.
    % tab delimiters between fields, comma delimiters between vector elements,
    % semicolons between matrix rows.
    
    %% concatenate immunity status variables to agent lists:
    % note: for hypotheticals this is disabled. 
    % to prepare for julia read-in
    
    %ensure datetime is imported correctly: 
    %opts = detectImportOptions(fname_residents_imstat);
    %opts = setvaropts(opts,'outbreak_date','InputFormat','dd/MM/yyyy');
    
    %res_imstat = readtable(fname_residents_imstat, opts);
    %res_imstat = table2struct(res_imstat);
    
    %ensure datetime is imported correctly: 
    %opts = detectImportOptions(fname_workers_imstat);
    %opts = setvaropts(opts,'outbreak_date','InputFormat','dd/MM/yyyy');
    %workers_imstat = readtable(fname_workers_imstat, opts);
    
    %n_med = size(workers_M, 2);
    %n_gen = size(workers_G, 2);
    %workers_m_imstat = table2struct(workers_imstat(1:n_med, :)); 
    %workers_g_imstat = table2struct(workers_imstat(n_med+1:end, :));
    
    % exclude staff_id and outbreak_index fields from imstat 
    % this avoids duplicate fields upon aggregation of properties
    %workers_m_imstat = rmfield(workers_m_imstat, 'staff_id');
    %workers_g_imstat = rmfield(workers_g_imstat, 'staff_id');
    %res_imstat = rmfield(res_imstat, 'resident_id');
    %workers_m_imstat = rmfield(workers_m_imstat, 'outbreak_index');
    %workers_g_imstat = rmfield(workers_g_imstat, 'outbreak_index');
    %res_imstat = rmfield(res_imstat, 'outbreak_index');
    
    workers_G_full = workers_G;
    workers_M_full = workers_M;
    residents_full = residents;
    
    
    %%
    write_agents_to_csv_Julia(residents_full, [outdir_fac_i '\residents_for_Julia_test']);
    write_agents_to_csv_Julia(workers_G_full, [outdir_fac_i '\workers_G_for_Julia_test']);
    write_agents_to_csv_Julia(workers_M_full, [outdir_fac_i '\workers_M_for_Julia_test']);
    

    % disabling all network generators - these are now constructed in
    % Julia from the agent-room assignments. 
    % write neighbour lists:
%     N_lists_struct = cell2struct(NeighbourLists, {'source_id', 'out_edges'}, 2);
%     N_lists_label = [outdir_fac_i '\Neighbour_lists_test' ];
%     write_network_to_csv_Julia(N_lists_struct, N_lists_label)
%     
%     %write weighted edge lists
%     E_list_struct = cell2struct(num2cell(Weighted_Networks), {'source_id', 'target_id', 'day', 'weight'}, 2);
%     E_list_label = [outdir_fac_i '\Edge_lists_test'];
%     write_network_to_csv_Julia(E_list_struct, E_list_label);
%     
%     %write multigraph edge lists
%     MGraph_struct = cell2struct(num2cell(Multigraph_EL), {'source_id', 'target_id', 'day', 'weight'}, 2);
%     MGraph_label = [outdir_fac_i '\MGraph_test'];
%     write_network_to_csv_Julia(MGraph_struct, MGraph_label);
    
    
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




