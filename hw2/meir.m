x =[1,0]'; %put x = [1;0] for (a), x = [-1;1] for (b)


%initial values of f(x_k), J(x_k)
fx = [(x(1))^2 - (x(2))^2 ; 2*(x(1))*(x(2)) - 1]; %f(x) where x is the initial vector
J = [2*(x(1)), -2*(x(2)); 2*(x(2)), 2*(x(1))]; %(J_f)(x) where x is the initial vector

k=0;

determinant = J(1,1) * J(2,2) - J(1,2) * J(2,1);
Data = [x(1,1), x(2,1), norm(fx), J(1,1), J(1,2), J(2,1), J(2,2)];
% [k, x', f',norm(f',2),J(1,1:2);'','','','','','','',J(2,1:2)];

%The matrix "Data" stores the index of iteration (k), x_k (as a row vector),
%f(x_k) (as a row vector), the 2 norm of f(x_k), the determinant of (J_f)(x_k), 
%and the first and second row of J(x_k), respectively given that we put the
%second row of J(x_k) by itself on a new line.

%NMS is a function I wrote to perform Newton's method for systems (specifically for this problem)


while (norm(fx,2) >= 10^-10) && (k < 30) && (det(J) ~= 0) 
    [x1, x2, fx,Jf] = NewtonM(x(1), x(2));
    k=k+1;
    
    n=norm(fx);
    determinant = Jf(1,1) * Jf(2,2) - Jf(1,2) * Jf(2,1);
    x = [x1, x2]';
    Data = [Data; x1, x2, n, Jf(1,1), Jf(1,2), Jf(2,1), Jf(2,2)];
%         k, x', f',norm(f',2), J(1,1:2);'','','','','','','',J(2,1:2)];
    %line 24 gives the logic to update "Data" with the data from the last
    %iteration of Newton's Method
end




