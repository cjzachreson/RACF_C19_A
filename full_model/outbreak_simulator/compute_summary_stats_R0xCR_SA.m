% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023
% updated on June 11th, 2024



% analysis of summary statistics from outbreak simulations


%% set up parameter sweep

clear all
close all

root_dir = pwd();
        
data_dir = fullfile(root_dir, 'output_v11_OB_SA_R0', 'facID_10000003_hyp');

output_dir = fullfile(root_dir, 'analysis', 'output_summary_stats_R0xCR_SA');

if ~isfolder(output_dir)
    mkdir(output_dir)
end

testing_strategies = ["asymp_testing_OB_only",...
                       ];

R0 = ["1.0", "1.5", "2.0", "2.4"];

bkgCR = ["3.0", "6.0", "10.0", "20.0", "30.0"];
     

R0xCR_pairs =...
    ["(1.0, 3.0)",... 
    "(1.0, 6.0)",... 
    "(1.0, 10.0)",... 
    "(1.0, 20.0)",...
    "(1.0, 30.0)",... 
    "(1.5, 3.0)",... 
    "(1.5, 6.0)",... 
    "(1.5, 10.0)",... 
    "(1.5, 20.0)",...
    "(1.5, 30.0)",... 
    "(2.0, 3.0)",... 
    "(2.0, 6.0)",... 
    "(2.0, 10.0)",... 
    "(2.0, 20.0)",... 
    "(2.0, 30.0)",... 
    "(2.4, 3.0)",... 
    "(2.4, 6.0)",... 
    "(2.4, 10.0)",... 
    "(2.4, 20.0)",...
    "(2.4, 30.0)"] ;

tscale_vals = [...
    "0.065",...
    "0.055",...
    "0.04",...
    "0.025",...
    "0.02",...
    "0.105",...
    "0.08",...
    "0.065",...
    "0.045",...
    "0.03",...
    "0.155",...
    "0.12",...
    "0.09",...
    "0.06",...
    "0.045",...
    "0.2",...
    "0.15",...
    "0.11",...
    "0.075",...
    "0.055"  ];

R0xCR_to_tscale = containers.Map(R0xCR_pairs, tscale_vals);


           
lockdown_compliance_OB = ["0.9", "0.0"] 

%% set up output table for summary stats: 

%medians and 90th quantiles: 

% cumulative incidence
% outbreak duration
% peak FTE deficit
% cumulative detected resident cases
% time to outbreak declaration
% time to first case detection

VarNames_q50 = {
            'p_outbreak',...
            'cum_I_q50',...
            'OB_duration_q50',...
            'peak_FTE_def_q50',...
            'cum_Res_Iso_q50',...
            't_declare_q50',...
            't_first_case_q50'};
        
VarNames_q90 = {
            'p_outbreak',...
            'cum_I_q90',...
            'OB_duration_q90',...
            'peak_FTE_def_q90',...
            'cum_Res_Iso_q90',...
            't_declare_q90',...
            't_first_case_q90'};
        
VarTypes = {'double', 'double','double','double','double','double','double'};

%set up row names corresponding to each scenario: 
[A,B] = meshgrid(R0,bkgCR);
c = cat(2,A,B);
d = reshape(c,[],2);
d = strcat('R0_', d(:, 1), '_CR_' ,d(:, 2));

[C, D] = meshgrid(d, lockdown_compliance_OB);
e = cat(2, C, D);
f = reshape(e, [], 2);
f = strcat(f(:, 1), '_LD_', f(:, 2));

row_names = cellstr(f);

summary_stats_q50 = table('Size', [size(row_names, 1), size(VarNames_q50, 2)],...
                          'VariableTypes', VarTypes,...
                          'VariableNames', VarNames_q50,...
                          'RowNames', row_names) ;

summary_stats_q90 = table('Size', [size(row_names, 1), size(VarNames_q90, 2)],...
                          'VariableTypes', VarTypes,...
                          'VariableNames', VarNames_q90,...
                          'RowNames', row_names) ;

             
       
%% plots of total cumulative incidence. 
for i = 1:size(testing_strategies, 2)
    
    strat_i = char(testing_strategies(i));

    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end


    
    data_dir_i = fullfile(strat_i);

   
    for r_j = 1:size(R0)

        R0_j = R0(r_j)

        data_dir_ij = fullfile(data_dir_i, ['R0_', strrep(char(R0_j), '.', 'p')])

        for cr_k = 1:size(bkgCR)
    
            CR_k = bkgCR(cr_k);

            data_dir_ijk = fullfile(data_dir_ij, ['bkgCR_', strrep(char(CR_k), '.', 'p')]);


            R0xCR_pair = strcat("(", R0_j, ", ", CR_k, ")");

            tscale = R0xCR_to_tscale(R0xCR_pair)

            data_dir_ijk = fullfile(data_dir_ijk, ['tscale_', strrep(char(tscale), '.', 'p')]);


            for ld_l = 1:size(lockdown_compliance, 2)

                LD_l = lockdown_compliance(ld_l);
                
                data_dir_ijkl = fullfile(data_dir_ijk, ['lockdown_compliance_', strrep(char(LD_l), '.', 'p')] );


                full_data_dir = fullfile(data_dir, data_dir_ijkl);

                output_linelist_fname = fullfile(full_data_dir, 'output_linelist.csv');

                output_linelist = readtable(output_linelist_fname);
                
                is_outbreak = strcmp(output_linelist.OB_declared, 'true');
                
                I_tot_All = output_linelist.I_tot_res +...
                           output_linelist.I_tot_staff;
                
                I_tot_OB = output_linelist.I_tot_res(is_outbreak) +...
                           output_linelist.I_tot_staff(is_outbreak);
                       
              
                if strcmp(strat_i, 'unmitigated')
                   I_tot_OB = I_tot_All(I_tot_All > 80); %HACK: this just plots the big outbreaks (for comparison)
                   I_tot = I_tot_All;
                else
                    I_tot = I_tot_OB;
                    
                end
                
                p_outbreak = numel(I_tot_OB) / numel(I_tot_All);
                
                cum_I_q50 = quantile(I_tot_OB, 0.5);
                cum_I_q90 = quantile(I_tot_OB, 0.9);
                
                scenario_label = strcat('R0_', R0_j, '_CR_', CR_k, '_LD_', LD_l);
                
                
                summary_stats_q50([scenario_label], ["p_outbreak"]) = {p_outbreak};
                summary_stats_q90([scenario_label], ["p_outbreak"]) = {p_outbreak};
                
                summary_stats_q50([scenario_label], ["cum_I_q50"]) = {cum_I_q50};
                summary_stats_q90([scenario_label], ["cum_I_q90"]) = {cum_I_q90};
                
                figure(i)
                d = 5;
                h = histogram(I_tot, 'BinEdges', [0:d:210], 'normalization', 'probability');
                xlabel('cumulative cases')
                hold on
                
                I_cum = (h.BinEdges(1:end-1))';
                frequency = (h.Values)';
                
                hist_table_I_tot = table(I_cum, frequency);
                
                hist_fname = ['I_tot_hist_',...
                    char(strrep(scenario_label, '.', 'p')), '.csv'];
                
                writetable(hist_table_I_tot, [output_dir, hist_fname])
                
                
                                            
            end

        end
    end
    
    Legend = cell(1, 1);
    
    for j = 1:size(lockdown_compliance, 2)
        LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
        Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
    end
    legend(Legend) 
end



%% plots of outbreak duration: 
for i = 1:size(testing_strategies, 2)
    
    strat_i = char(testing_strategies(i));
    
    data_dir_i = [data_dir, strat_i, '\'];
    
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    
    
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        
        t_OB_on = output_linelist_ij.t_OB_on(is_outbreak);   
        
        t_OB_over = output_linelist_ij.sim_duration(is_outbreak);
        
        OB_duration = t_OB_over - t_OB_on;
               
        % output summary statistics. 
        OB_duration_q50 = quantile(OB_duration, 0.5);
        OB_duration_q90 = quantile(OB_duration, 0.9);
        
        scenario_label_ij = ...
            strcat(testing_strategies(i), '_LD_', lockdown_compliance(j));
        
        summary_stats_q50([scenario_label_ij], ["OB_duration_q50"]) = {OB_duration_q50};
        summary_stats_q90([scenario_label_ij], ["OB_duration_q90"]) = {OB_duration_q90};
        
        figure(i + size(testing_strategies, 2)*2)
        d = 5;
        h = histogram(OB_duration, 'BinEdges', [0:d:100], 'normalization', 'probability');
        xlabel('outbreak duration')
        hold on
        
        duration = (h.BinEdges(1:end-1))';
        frequency = (h.Values)';
        
        hist_table_duration = table(duration, frequency);
        
        hist_fname = ['OB_duration_hist_',...
                      char(strrep(scenario_label_ij, '.', 'p')), '.csv'];
        
        writetable(hist_table_duration, [output_dir, hist_fname])
        
                                    
    end
    
    Legend = cell(1, 1);
    
    if ~isempty(lockdown_compliance)
        for j = 1:size(lockdown_compliance, 2)
            LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
            Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
        end
        legend(Legend) 
    end
end



%% plot outbreak duration vs. cumulative cases
%{
for i = 1:size(testing_strategies, 2)
    
    strat_i = char(testing_strategies(i));
    
    data_dir_i = [data_dir, strat_i, '\'];
    
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    
    
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        
        t_OB_on = output_linelist_ij.t_OB_on(is_outbreak);   
        
        t_OB_over = output_linelist_ij.sim_duration(is_outbreak);
        
        OB_duration = t_OB_over - t_OB_on;
        
        
        I_tot_OB = output_linelist_ij.I_tot_res(is_outbreak) +...
                   output_linelist_ij.I_tot_staff(is_outbreak);
               
        figure(i + size(testing_strategies, 2)*2)
        scatter(I_tot_OB, OB_duration, '.')
        ylabel('outbreak duration vs. cumulative incidence')
        hold on
        
                                    
    end
    
    Legend = cell(2, 1);
    
    for j = 1:size(lockdown_compliance, 2)
        LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
        Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
    end
    legend(Legend) 
end
%}



%% plot peak FTE deficit
for i = 1:size(testing_strategies, 2)
    strat_i = char(testing_strategies(i));
    data_dir_i = [data_dir, strat_i, '\'];
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        peak_FTE_def = output_linelist_ij.FTE_def_max(is_outbreak);   
        
        
        % output summary statistics. 
        peak_FTE_def_q50 = quantile(peak_FTE_def, 0.5);
        peak_FTE_def_q90 = quantile(peak_FTE_def, 0.9);
        
        scenario_label_ij = ...
            strcat(testing_strategies(i), '_LD_', lockdown_compliance(j));
        
        summary_stats_q50([scenario_label_ij], ["peak_FTE_def_q50"]) = {peak_FTE_def_q50};
        summary_stats_q90([scenario_label_ij], ["peak_FTE_def_q90"]) = {peak_FTE_def_q90};
        
        figure(i + size(testing_strategies, 2)*3)
        d = 2;
        h = histogram(peak_FTE_def, 'BinEdges', [0:d:50], 'normalization', 'probability');
        xlabel('peak FTE deficit')
        hold on
        
        max_FTE_def = (h.BinEdges(1:end-1))';
        frequency = (h.Values)';
        
        hist_table_FTE = table(max_FTE_def, frequency);
        
        hist_fname = ['FTE_def_hist_',...
                      char(strrep(scenario_label_ij, '.', 'p')), '.csv'];
        
        writetable(hist_table_FTE, [output_dir, hist_fname])
        
                                    
    end
    
    if ~isempty(lockdown_compliance)
        Legend = cell(1, 1);

        for j = 1:size(lockdown_compliance, 2)
            LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
            Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
        end
        legend(Legend) 
    end
end



%% plot cumulative residents detected (isolated)
for i = 1:size(testing_strategies, 2)
    strat_i = char(testing_strategies(i));
    data_dir_i = [data_dir, strat_i, '\'];
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        cum_Res_Iso = output_linelist_ij.Det_tot_res(is_outbreak);   
        
        
        % output summary statistics. 
        cum_Res_Iso_q50 = quantile(cum_Res_Iso, 0.5);
        cum_Res_Iso_q90 = quantile(cum_Res_Iso, 0.9);
        
        scenario_label_ij = ...
            strcat(testing_strategies(i), '_LD_', lockdown_compliance(j));
        
        summary_stats_q50([scenario_label_ij], ["cum_Res_Iso_q50"]) = {cum_Res_Iso_q50};
        summary_stats_q90([scenario_label_ij], ["cum_Res_Iso_q90"]) = {cum_Res_Iso_q90};
        
        figure(i + size(testing_strategies, 2)*4)
        d = 5;
        h = histogram(cum_Res_Iso, 'BinEdges', [0:d:90],...
            'normalization', 'probability');
        xlabel('total isolated residents (detected)')
        hold on
        
        tot_Res_Iso = (h.BinEdges(1:end-1))';
        frequency = (h.Values)';
        
        hist_table_Iso = table(tot_Res_Iso, frequency);
        
        hist_fname = ['Iso_hist_',...
                      char(strrep(scenario_label_ij, '.', 'p')), '.csv'];
        
        writetable(hist_table_Iso, [output_dir, hist_fname])
        
                                    
    end
    
    if ~isempty(lockdown_compliance)
        Legend = cell(1, 1);

        for j = 1:size(lockdown_compliance, 2)
            LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
            Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
        end
        legend(Legend) 
    end
end



%% plot time to outbreak declaration
for i = 1:size(testing_strategies, 2)
    strat_i = char(testing_strategies(i));
    data_dir_i = [data_dir, strat_i, '\'];
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        t_outbreak_declared = output_linelist_ij.t_OB_on(is_outbreak);   
        
        
        % output summary statistics. 
        t_declare_q50 = quantile(t_outbreak_declared, 0.5);
        t_declare_q90 = quantile(t_outbreak_declared, 0.9);
        
        scenario_label_ij = ...
            strcat(testing_strategies(i), '_LD_', lockdown_compliance(j));
        
        summary_stats_q50([scenario_label_ij], ["t_declare_q50"]) = {t_declare_q50};
        summary_stats_q90([scenario_label_ij], ["t_declare_q90"]) = {t_declare_q90};
        
        figure(i + size(testing_strategies, 2)*5)
        d = 2;
        h = histogram(t_outbreak_declared, 'BinEdges', [0:d:90],...
            'normalization', 'probability');
        xlabel('delay until outbreak declaration')
        hold on
        
        t_OB_declare = (h.BinEdges(1:end-1))';
        frequency = (h.Values)';
        
        hist_table_declare = table(t_OB_declare, frequency);
        
        hist_fname = ['t_declare_hist_',...
                      char(strrep(scenario_label_ij, '.', 'p')), '.csv'];
        
        writetable(hist_table_declare, [output_dir, hist_fname])
        
                                    
    end
    
    if ~isempty(lockdown_compliance)
        Legend = cell(1, 1);

        for j = 1:size(lockdown_compliance, 2)
            LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
            Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
        end
        legend(Legend) 
    end
end



%% plot time to first detected case
for i = 1:size(testing_strategies, 2)
    strat_i = char(testing_strategies(i));
    data_dir_i = [data_dir, strat_i, '\'];
    if strcmp(strat_i, 'unmitigated')
        lockdown_compliance = [];
    else
        lockdown_compliance = lockdown_compliance_OB;
    end
    for j = 1:size(lockdown_compliance, 2)
       
        LD_j = strrep(char(lockdown_compliance(j)), '.', 'p');
        
        data_dir_ij = [data_dir_i, 'lockdown_compliance_' LD_j '\'];
        
        output_linelist_ij = readtable([data_dir_ij,...
                                        'output_linelist.csv']);
        
        is_outbreak = strcmp(output_linelist_ij.OB_declared, 'true');
        
        t_first_detection = output_linelist_ij.t_first_detection(is_outbreak);   
        
        
        % output summary statistics. 
        t_first_case_q50 = quantile(t_first_detection, 0.5);
        t_first_case_q90 = quantile(t_first_detection, 0.9);
        
        scenario_label_ij = ...
            strcat(testing_strategies(i), '_LD_', lockdown_compliance(j));
        
        summary_stats_q50([scenario_label_ij], ["t_first_case_q50"]) = {t_first_case_q50};
        summary_stats_q90([scenario_label_ij], ["t_first_case_q90"]) = {t_first_case_q90};
        
        figure(i + size(testing_strategies, 2)*6)
        d = 2;
        h = histogram(t_first_detection, 'BinEdges', [0:d:90],...
            'normalization', 'probability');
        xlabel('delay until first detected case')
        hold on
        
        t_first_case = (h.BinEdges(1:end-1))';
        frequency = (h.Values)';
        
        hist_table_detect = table(t_first_case, frequency);
        
        hist_fname = ['t_declare_hist_',...
                      char(strrep(scenario_label_ij, '.', 'p')), '.csv'];
        
        writetable(hist_table_detect, [output_dir, hist_fname])
        
                                    
    end
    
    if ~isempty(lockdown_compliance)
        Legend = cell(1, 1);

        for j = 1:size(lockdown_compliance, 2)
            LD_j = strrep(num2str(lockdown_compliance(j)), '.', 'p');
            Legend{j} = [strrep(strat_i, '_', ' ') ',  LD: ' strrep(LD_j, 'p', '.')];
        end
        legend(Legend) 
    end
end

%finally, print table of quantiles: 

writetable(summary_stats_q50, [output_dir, 'summary_stats_q50.csv'],'WriteRowNames',true)   
writetable(summary_stats_q90, [output_dir, 'summary_stats_q90.csv'],'WriteRowNames',true)   
