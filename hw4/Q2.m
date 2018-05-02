x = [10, 1];
t = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
y = [6.3, 3.3, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15];
a = 1;

yNew = log(y);
A = ones(2,8);
A(2,:) = t;
A = A';

x = inv(A' * A) * A' * yNew;


for index = 1:1:100
	err_x1 = 0;
	err_x2 = 0;
	for i = 1:1:8
		err_x1 = err_x1 + 2* x(1) * (exp(x(2) * t(i)))^2 - 2*exp(x(2) * t(i))*y(i);
		err_x2 = err_x2 + 2* (x(1) ^ 2) * t(i) * (exp(x(2) * t(i)))^2 - 2*x(1) * t(i) * exp(x(2) * t(i))*y(i);

	end
	delta_fx = [err_x1, err_x2];
	x
	delta_fx
	x = x - a*delta_fx;
	a = a/2;
end

fx = 0;
for i = 1:1:8
	fx = fx + (x(1)*exp(x(2) * t(i)) - y(i))^2;
end

fx
&