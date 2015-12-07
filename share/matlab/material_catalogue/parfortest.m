function parfortest()
pause on
tic;
N=500;
for i=1:N
    sequential_answer=slow_fun(i);
end
sequential_time=toc
tic;
parfor i=1:N
   sequential_answer=slow_fun(i);
end
parallel_time=toc
end

function result=slow_fun(x)
    pause(0.01);
    result=x;
end