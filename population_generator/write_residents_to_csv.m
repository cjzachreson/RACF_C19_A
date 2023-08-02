function f = write_residents_to_csv(residents, label)

fname = [label, '.csv'];

output_string = [];

c_names = fieldnames(residents);

delim = {']['};

linestring_0 = ['[', strjoin(c_names, delim), ']'];

output_string = [output_string; linestring_0];

for i = 1:numel(residents)

linestring_i = [];

res_i = residents(i);

res_i_c = struct2cell(res_i);

for j = 1:size(res_i_c, 1)

        
    cell_elements = res_i_c{j};
    
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








