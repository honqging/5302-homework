% x1 = c(1), x2 = c(2)
function [hessian] = getErrHessian(y,t,c)
	hessian = [0,0;0,0];
	for i = 1:1:8
		hessian(1, 1) = hessian(1, 1) + 2 * (exp(c(2) * t(i))) ^ 2;
		hessian(1, 2) = hessian(1, 2) + 4 * c(1) * t(i) * exp(c(2) * t(i));
		hessian(2, 1) = hessian(1, 2);
		hessian(2, 2) = hessian(2, 2) + 4 * (c(1)) * 2 * (t(i)) ^ 2 * (exp(c(2) * t(i))) ^ 2 - 2 * c(1) * (t(i)) ^ 2 * y(i) * exp(c(2) * t(i));
		
	end
end