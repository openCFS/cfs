function res = hillmandeltovoigt(tensor)

Q = zeros(3,3);

Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1/sqrt(2);

res = Q' * tensor * Q;