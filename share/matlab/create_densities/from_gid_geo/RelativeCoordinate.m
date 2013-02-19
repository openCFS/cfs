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

elseif abs(rDeltaZ) > abs(rDeltaY) % && abs(rDeltaZ) > abs(rDeltaY)
    if rDeltaZ == 0
       return;
    end

    rLambda = (p(3) - s(3)) / rDeltaZ;
    return;

else
    if rDeltaY == 0
        return;
    end

    if abs(rDeltaX) == abs(rDeltaY) && abs(rDeltaX) == abs(rDeltaZ)
      d1 = p(1) - s(1);
      d2 = p(2) - s(2);
      d3 = p(3) - s(3);
      if abs(d1) > abs(d2) && abs(d1) > abs(d3)
        rLambda = d1 / rDeltaX;
        return;
      elseif abs(d3) > abs(d2)
        rLambda = d3 / rDeltaZ;
        return;
      else 
        rLambda = d2 / rDeltaY;
        return;
      end
      end
    end

    rLambda = (p(2) - s(2)) / rDeltaY;
    return;
end
end
