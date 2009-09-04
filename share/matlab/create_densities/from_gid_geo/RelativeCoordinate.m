function rLambda = RelativeCoordinate (s, e, p)

rDeltaX = e(1) - s(1);
rDeltaY = e(2) - s(2);
rDeltaZ = e(3) - s(3);
rLambda = 0;

if abs(rDeltaX) > abs(rDeltaY) && abs(rDeltaX) > abs(rDeltaZ)
    if rDeltaX == 0
        return;
    end

    rLambda = (p(1) - s(1)) / rDeltaX;
    return;

elseif abs(rDeltaZ) > abs(rDeltaY) && abs(rDeltaZ) > abs(rDeltaY)
    if rDeltaZ == 0
       return;
    end

    rLambda = (p(3) - s(3)) / rDeltaZ;
    return;

else
    if rDeltaY == 0
        return;
    end

    rLambda = (p(2) - s(2)) / rDeltaY;
    return;
end
end