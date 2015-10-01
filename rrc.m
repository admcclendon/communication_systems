function [ o ] = rrc( B, Ts, t )
%RRC Summary of this function goes here
%   Detailed explanation goes here

    if (t == 0)
        o = 1 - B + 4*(B/pi);
    elseif (t == Ts/(4*B) || t == -Ts/(4*B))
        o = (B/sqrt(2))*((1+2/pi)*sin(pi/(4*B))+(1-2/pi)*cos(pi/(4*B)));
    else
        o = (sin(pi*(t/Ts)*(1-B))+4*B*(t/Ts)*cos(pi*(t/Ts)*(1+B)))/(pi*(t/Ts)*(1-(4*B*(t/Ts))^2));
    end
end

