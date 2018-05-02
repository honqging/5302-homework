t = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
y = [6.3, 3.3, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15];

c = [10; -1];
for idx = 1:1:100
	grad = getErrGradient(y, t, c);
	hessian = getErrHessian(y, t, c);
	s0 = hessian' * (-grad);
	c = c + s0;
end



