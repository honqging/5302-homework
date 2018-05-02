k = [2 4 6 8 10];

% k = 2, error = 8.9685e-16
maxError = [0.0000 0.0331 0.0528 0.0684 0.0811];

plot(k,maxError,'r-o');
title('relation between maximum absolute differences and k');
xlabel('k');
ylabel('maximum absolute difference')


