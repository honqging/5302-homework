% x1 = c(1), x2 = c(2)
function [gradient] = getErrGradient(y,t,c)
	gradient = [0;0];
	for i = 1:1:8
		gradient(1) = gradient(1) + 2* c(1) * (exp(c(2) * t(i))) ^ 2 - 2 * exp(c(2) * t(i));
		gradient(2) = gradient(2) + 2 * (c(1))^2 * t(i) * (exp(c(2) * t(i))) ^ 2 - 2 * c(1) * t(1) * exp(c(2) * t(i)) * y(i);
	end
end 