[x1, x2, Jf] = NewtonM(x1, x2);


fx = [x1 ^ 2 - x2 ^ 2; 2 * x1 * x2 - 1];
n=norm(fx);
determinant = Jf(1,1) * Jf(2,2) - Jf(1,2) * Jf(2,1);

display([x1, x2]);
display(Jf);
display(n);
display(determinant)
        




