% Author: Cameron Zachreson
% Institution: The University of Melbourne
% Simulation code acompanying the manuscript entitled: 
% "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
% Date released: Dec. 18, 2023

function numeric_output = workers_from_csv(workers_fname)


fid = fopen(workers_fname, 'rt');

field_names = [];

delim = {'][', '] [', ']\n ['};

test = textscan(fid, '%s', 'Delimiter', delim) ;

fclose(fid);

sz = size(test{1,1}, 1);
% remove dangling brackets
for i = 1:sz
    test{1,1}{i} = strrep(test{1, 1}{i}, ']', '');
    test{1,1}{i} = strrep(test{1, 1}{i}, '[', '');
end



%determine the number of columns by reading field names
n_col = 0;
for i = 1:sz
    if ~isnan(str2double(test{1,1}{i}))
        break
    end
    n_col = n_col + 1;
    
end


% reshape
n_row = sz/n_col;
test_2d = {};
j = 1;
k = 1;
for i = 1:sz
    
    test_2d{k, j} = test{1, 1}{i};
    
    j = j + 1;
    if j == n_col + 1
        j = 1;
        k = k + 1;
    end
end

field_names = {test_2d{1, :}};

% rooms_text = struct();
%
% for i = 1:n_col
%
%     field_name = field_names{i};
%
%     for j = 1:n_row-1
%
%     field_vals = test_2d{j+s1, i};
%
%     rooms_text(j).(field_name) = field_vals;
%
%     end
% end

numeric_output = struct;
%convert to numeric:
for i = 1:n_col
    
    field_name = field_names{i};
    
    for j = 1:n_row-1
        
        if strcmp(field_name, 'outbreak_key')
            field_vals = test_2d{j+1, i};
        else
            
            % use of matlab delimiters in text file comes in handy here.
            field_vals = eval([ '[' test_2d{j+1, i} ']' ]);
            
        end
        
        numeric_output(j).(field_name) = field_vals;
        
    end
    
end











