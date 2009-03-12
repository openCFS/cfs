# extract the change of design 
# input:  Array with "id data1 data2 data3 ..." as lines
# output: Array with the id as first column 
#         and normalized design change to prior as second column.
#         Output misses a line with the first id as there is no prior data
function retval = opt_design_change (A)
  [rows,cols] = size(A);
  # our result is one line shorter
  retval = zeros(rows-1, 2);
  for current = 2:rows
    last = current-1;
    # our result is one line shorter, such the index is last 
    retval(last,1) = A(current, 1);
    # sum up the differences
    diff = 0;
    for inner = 2:cols
      diff += abs(A(current, inner) - A(last, inner));
    endfor     
    diff /= (cols-1);
    retval(last,2) = diff;
  endfor
endfunction