function hom_stats_to_mat(inputfile,outputfile,s1,s2)

%read material catalogue
list = load(inputfile);
A = zeros(3,3);
for it = 1:length(s1)
    for jt = 1:length(s2)
        for i=2:size(list,1)
            if s1(it) == list(i,1) && s2(jt) == list(i,2) 
                A(1,1) = list(i,3);
                A(1,2) = list(i,4);
                A(2,1) = list(i,4);
                A(2,2) = list(i,5);
                A(3,3) = list(i,6);
            end
        end
        %dlmwrite('A.mat',A,' ',' ','1e-9','%5.8g');
    end
end
dlmwrite(outputfile, A, 'delimiter', ',', 'precision', 9); 
