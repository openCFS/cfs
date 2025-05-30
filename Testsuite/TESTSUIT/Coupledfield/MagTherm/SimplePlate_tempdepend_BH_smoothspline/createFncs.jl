#%%
using QuadGK
using DelimitedFiles

f100(B_R) = 39.2004*exp(0.00094279*B_R^(11.5695))-21.1455
f110(B_R) = 0.0001666*exp(9.7377*B_R^(0.88664))+29.3703
f120(B_R) = 1.2733e-05*exp(13.9035*B_R^(0.7816))+57.7738
f130(B_R) = 9.8337e-06*exp(15.168*B_R^(0.82632))+122.7577
f140(B_R) = 0.00031078*exp(11.5428*B_R^(1.1784))+195.7396

function evalInt(f, i)
    return quadgk(f, 0, i, order=13)[1]
end

function write2Fnc(f, B, filename)
    data = [[evalInt(f, i) for i in B] collect(B)]
    writedlm(filename, data, ' ')
end


B = range(0.0 ,2.0, 15)


write2Fnc(f100, B, "BH100.fnc")
write2Fnc(f110, B, "BH110.fnc")
write2Fnc(f120, B, "BH120.fnc")
write2Fnc(f130, B, "BH130.fnc")
write2Fnc(f140, B, "BH140.fnc")
