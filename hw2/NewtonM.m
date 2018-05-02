function [y1, y2, fx, Jf] = NewtonM(x1, x2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    fx = [x1 ^ 2 - x2 ^ 2; 2 * x1 * x2 - 1];
    Jf = [2 * x1, -2 * x2; 2 * x2, 2 * x1];
    
    s0 = Jf\-fx;
    
    y1 = x1 + s0(1,1);
    y2 = x2 + s0(2,1);
end

