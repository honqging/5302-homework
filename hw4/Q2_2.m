t = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]';
y = [6.3, 3.3, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';

T = (sqrt(5)-1)/2;
a = -1.38;
b = -1.2;

c = [a + (1-T)*(b-a), a + T*(b-a)];

err = zeros(1,2);


for j=1:2
	
	for i = 1:1:8
		A = exp(c(j) * t);
		x1 = inv(A' * A) * A' * y;
		err(:,j) = err(:,j) + (x1*exp(c(j) * t(i)) - y(i))^2;
	end
end

while b-a > 10^(-7)

	if(err(1) > err(2))
		a = c(1);
		c(1) = c(2);
		err(1) = err(2);
		c(2) = a + T*(b-a);
		err(2) = getErr(y, t, c(2));
	else
		b = c(2);
		c(2) = c(1);
		err(2) = err(1);
		err(1) = getErr(y,t,c(1));
	end 
end


