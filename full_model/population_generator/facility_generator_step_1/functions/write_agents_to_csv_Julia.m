% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function write_agents_to_csv_Julia(agents, label)

fname = [label, '.csv'];

output_string = [];

c_names = fieldnames(agents);

delim = {'\t'};

linestring_0 = strjoin(c_names, delim);

output_string = [output_string; linestring_0];

for i = 1:numel(agents)

linestring_i = [];

agent_i = agents(i);

agent_i_c = struct2cell(agent_i);

for j = 1:size(agent_i_c, 1)

        
    cell_elements = agent_i_c{j};
    
    line_entry_in = string(cell_elements);
    
    if numel(line_entry_in) > 1 && isvector(line_entry_in) %vector-valued properties (e.g., roster)
    
        line_entry_out = strjoin(line_entry_in, ',');
        
    elseif numel(line_entry_in) > 1 && ismatrix(line_entry_in) %vector-valued properties (e.g., roster)
        
        num_rows = size(line_entry_in, 1);
        
        line_entry_out = [];
        
        for k = 1:num_rows
            
            if isempty(line_entry_out)
                line_entry_out = strjoin(line_entry_in(k, :), ',');
            else
            
            line_entry_out = [line_entry_out, [';', strjoin(line_entry_in(k, :), ',')] ];
            end
        end
        
        
        line_entry_out = strjoin(line_entry_out);
        
    else %single value
        line_entry_out = line_entry_in;
        
    end
    
    if isempty(linestring_i)
        linestring_i = line_entry_out;
    else
    
    linestring_i = [linestring_i, [delim{1}, line_entry_out]];
    
    end
end

output_string = [output_string; strjoin(linestring_i, '')];

end

output_string = cellstr(output_string);
output_string = strjoin(output_string', '\n');

fid = fopen(fname, 'wt');
fprintf(fid, output_string);
fclose(fid);








