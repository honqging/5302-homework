function y = water(a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    x = 11;

    g = 9.8065;
    r = 0.3;
    v = 15;

    y = x * (g/(r * v * cosd(a)) + tand(a)) + g/(r^2) * log(1 - r * x/(v * cosd(a)));

end

