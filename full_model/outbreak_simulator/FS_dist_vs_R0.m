% plots final size dist for each kappa, R0, for presentation as an image

input_data_dir = ['C:\Users\czachreson\Desktop\compositions_in_progress\'...
                   'RACF_A_outbreak_response\code\RACF_C19_A_repo\full_model',...
                   '\outbreak_simulator\output_v9_FS_test_3\immunity_off\',...
                   'facID_10000003_hyp\n_100'];
               
               
kappa = 0.005:0.005:0.5;

xmax = 210;

bin_width = 2;

bin_edges = [0:bin_width:xmax];

output_mat = NaN(numel(bin_edges)-1, numel(kappa));

for i = 1:numel(kappa)
    
    k = kappa(i);
    
    input_fname = [input_data_dir, '\tscale_' strrep(num2str(k), '.', 'p' ),...
                  '\output_summary.csv'];
              
    FS = readtable(input_fname, 'ReadVariableNames', true);
    FS = FS.secondary_cases_tot;
    
    FS_hist = histogram(FS, 'BinEdges', bin_edges, 'normalization', 'probability');
    
    p_FS = FS_hist.Values;
    
    output_mat(:, i) = p_FS;
    
end

dlmwrite('analysis\p_FS_vs_kappa_mat.csv', output_mat')
              
    