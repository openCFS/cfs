function [x] = createFrame(nelx,nely,a,b)

x = ones(nelx,nely);

a = a/2;
b = b/2;

nela = round(a*nelx);
nelb = round(b*nely);

for elx = nela+1:(nelx-nela)
    for ely = nelb+1:(nely-nelb)
        x(elx,ely) = 1e-5;
    end
end
