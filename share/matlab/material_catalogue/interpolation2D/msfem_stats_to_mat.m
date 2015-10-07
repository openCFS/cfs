function msfem_stats_to_mat(inputfile,outputfile,s1,s2)

%read material catalogue
list = load(inputfile);
A = zeros(8,8);
row = 1;
column = 1;
count = 1.;
for it = 1:length(s1)
    for jt = 1:length(s2)
        for i=2:size(list,1)
            if s1(it) == list(i,1) && s2(jt) == list(i,2) 
                for j=3:38
                    A(row,column) = list(i,j);
                    A(column,row) = list(i,j);
                    column = column + 1;
                    if column > 8
                        count = count +1;
                        column = count;
                        row = row + 1;
                    end
                end
            end
        end
        %dlmwrite('A.mat',A,' ',' ','1e-9','%5.8g');
    end
end
dlmwrite(outputfile, A, 'delimiter', ',', 'precision', 9); 
