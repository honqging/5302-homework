alpha = zeros(7, 12);
b= zeros(7,12);

b(1,:) = [4 3 1 0 -1 -3 -4 -3 -1 0 1 3];
b(2,:) = [3 4 3 1 0 -1 -3 -4 -3 -1 0 1];
b(3,:) = [1 3 4 3 1 0 -1 -3 -4 -3 -1 0];
b(4,:) = [0 1 3 4 3 1 0 -1 -3 -4 -3 -1];
b(5,:) = [3 3 1 1 -1 -1 -3 -3 -1 -1 1 1];
b(6,:) = [3 1 -1 -3 -1 1 3 1 -1 -3 -1 1];
b(7,:) = [2 0 -2 0 2 0 -2 0 2 0 -2 0];


i = sqrt(-1);
n = 12;

A = zeros(12, 12);

for x = 1 : 12
	for k = 1 : 12
		A(x, k) = 1/sqrt(n) * exp(2 * pi * i * (k-1) * (x-1) / n);
	end
end

alpha = inv(A' * A) * A' * b';



