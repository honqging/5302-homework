function [err] = getErr(y,t,c)
	err = 0;
	for i = 1:1:8
		A = exp(c * t);
		x1 = inv(A' * A) * A' * y;
		err = err + (x1*exp(c * t(i)) - y(i))^2;
	end
end

