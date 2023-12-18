% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function f = write_workers_to_csv(workers, label)

fname = [label, '.csv'];

output_string = [];

c_names = fieldnames(workers);

delim = {']['};

linestring_0 = ['[', strjoin(c_names, delim), ']'];

output_string = [output_string; linestring_0];

for i = 1:numel(workers)

linestring_i = [];

w_i = workers(i);

w_i_c = struct2cell(w_i);

for j = 1:size(w_i_c, 1)

        
    cell_elements = w_i_c{j};
    
    line_entry_in = string(cell_elements);
    
    if numel(line_entry_in) > 1 && isvector(line_entry_in) %vector-valued properties (e.g., roster)
    
        line_entry_out = strjoin(['[' strjoin(line_entry_in, ',') ']'], '');
        
    elseif numel(line_entry_in) > 1 && ismatrix(line_entry_in) %vector-valued properties (e.g., roster)
        
        num_rows = size(line_entry_in, 1);
        
        line_entry_out = ['['];
        
        for k = 1:num_rows
            
            line_entry_out = [line_entry_out, strjoin([strjoin(line_entry_in(k, :), ','), ';'], '' ) ];
        
        end
        
        line_entry_out = [line_entry_out, ']'];
        
        line_entry_out = strjoin(line_entry_out);
        
    else %single value
        line_entry_out = strjoin(['[', line_entry_in, ']'], '');
        
    end
        
    linestring_i = [linestring_i, line_entry_out];
    
    
end

output_string = [output_string; strjoin(linestring_i, '')];

end

output_string = cellstr(output_string);
output_string = strjoin(output_string', '\n');

fid = fopen(fname, 'wt');
fprintf(fid, output_string);
fclose(fid);








