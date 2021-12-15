function [dEda] = monocubic_partialderiv(a, E, deriv_a)
% Approximation der Ableitungen in den Interpolationspunkten mit Finiten Differenzen
% vectors deriv_a and deriv_b are fixed derivatives used by the penalization
% approach
dEda = zeros(length(a),1);

% Central difference quotient for inner points
for j=2:length(a)-2
    dEda(j,:)=(E(j+1,:)-E(j-1,:))/(a(j+1)-a(j-1));
end

% Backwards difference quotient for outer points and next to outer points.
% As a=0 and a=1 are treated as full material, these points can
% not be used for central difference quotients.
dEda(1,:) = (E(2,:)-E(1,:))/(a(2)-a(1));
if length(a) > 2
    dEda(length(a)-1,:) = (E(length(a)-1,:)-E(length(a)-2,:))/(a(length(a)-1)-a(length(a)-2));
end

if isempty(deriv_a)
    dEda(length(a),:) = (E(length(a),:)-E(length(a)-1,:))/(a(length(a))-a(length(a)-1));
else
    dEda(length(a),:) = deriv_a(:);
end
