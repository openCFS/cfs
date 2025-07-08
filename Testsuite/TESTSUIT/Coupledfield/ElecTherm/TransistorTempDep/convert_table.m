function convert_table(fname,n_slice)

% Function to convert .table data from KAI to .dat files needed for the
% tripolemodel. The data is not scaled so we assume a factor of 1 for the
% definition within the material file (different to the one in the
% Testsuite example!). Values with up to one decimal place have been tested
% and will be written in the format XVX. Example: 1.5 V --> 1V5.
% Important: The function expects the specific resistance as an input
% parameter for the table conversion, since it inverts the value to get the
% conductivity in line 68.

% Input:
% fname ... file name of the .table file
% n_slice ... number of slices where VDS is constant

% Output: 
% Vds_const_slices.dat ... .dat file containing the list of slices
% slice_XXVX.dat ... data slice where VDS is constant. Each slice contains
% blocks where VGS is constant and the conductivity is given as a function
% of the temperature.

disp(strcat("Converting ", fname))

% Read data and init variables
data = readmatrix(fname);

[n,m] = size(data);

input_block_height = n/n_slice;

vds_vec = zeros(n_slice,1);
slice_names_vec = cell(n_slice,1);
vgs_vec = zeros(input_block_height-1,1);

% loop over all blocks (all values for vds)
for ii=1:n_slice
    % create a name for the current slice and get the Vds value
    vds_vec(ii,1) = data(input_block_height*(ii-1)+1,1);
    % prep the name
    vds_str = num2str(vds_vec(ii,1));
    if contains(vds_str,'.')
        % we have to split it
        vds_str_split = strsplit(vds_str,'.');
        cur_name = strcat('slice_',vds_str_split{1},'V',vds_str_split{2},'.dat');
        slice_names_vec{ii} = cur_name;
    else
        % easier version
        cur_name = strcat('slice_',vds_str,'V0','.dat');
        slice_names_vec{ii} = cur_name;
    end
    
    
    % loop over all values for vgs
    % create file
    fileID = fopen(slice_names_vec{ii},'w');
    for uu=1:input_block_height-1
        %write data within output-block and loop over vgs changes
        vgs_vec(uu) = data(input_block_height*(ii-1)+uu+1,1); % +1 offset because the first entry of each block is actually Vds
        
        % reset output block
        output_block_data_mat = zeros(m-1,3);
        
        output_block_data_mat(1:m-1,1) = vgs_vec(uu);
        % loop over the temperature
        for oo=1:m-1
            % write data for current output-block (only temp changes)
            output_block_data_mat(oo,2) = data(input_block_height*(ii-1)+1,1+oo);
            output_block_data_mat(oo,3) = 1/data(input_block_height*(ii-1)+uu+1,1+oo); % we have to invert the value since we need the conductivity but got the specific resistance
            
        end
        % write output block
        fprintf(fileID,'%10.4e %10.4e %10.4e\n',transpose(output_block_data_mat));
        % space between the output-blocks
        fprintf(fileID,'\n');
    end
    %fclose(fileID);
end

% finally, create the file combining all slices
fileID = fopen('Vds_const_slices.dat','w');
for ii=1:n_slice
    fprintf(fileID,'%10.4e %s\n',vds_vec(ii),slice_names_vec{ii});
end
fclose(fileID);

disp("Done!")

end